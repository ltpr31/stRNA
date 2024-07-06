
workdir=/path/to/work/
script=$workdir/script



####### seurat object construction

for i in mock Mo6hpi Mo12hpi Mo24hpi ; do
Rscript $script/spatial_obj_construction.r  \
--path  $workdir/data/${i}/  \
--png   $workdir/data/${i}/tissue.png  \
-o      $workdir/obj/   \
-p ${i}
done


##### construction ALBST with LOONG v1.0
##install LOONG v1.0
sudo dpkg -i LOONG_v1.0.deb
for i in mock Mo6hpi Mo12hpi Mo24hpi ; do
cd  $workdir/data/${i}
gunzip -c barcodes_pos.tsv.gz | awk -v OFS='\t' '{temp = $2; $2 = $3; $3 = temp; print}' > position.txt
done
##### a example is as follow
/usr/local/LOONG/bin/draw_a_loong --rootPath  $workdir/data/Mo6hpi  \
--txtName position.txt  \
--imgName tissue.png \
--outTxtName albst.txt \
--row_std 2 --outImg_height 80  \
--bgImg_height 1005 --bgImg_width 1000


#note:ALBST values were scaled to 0 to 2 by parameter "--row_std 2",
#and 1 was then subtracted from each ALBST value to limit ALBST from -1 to 1

######### add ALBST value
### "albst.txt" used in the paper were saved in the path "--albst $workdir/data/${i}/albst.txt"
for i in mock Mo6hpi Mo12hpi Mo24hpi ; do
Rscript $script/spatial_add_albst.r  \
--albst $workdir/data/${i}/out/albst.txt  \
-r $workdir/obj/${i}.rds  \
-o $workdir/obj_albst  \
-p ${i}  --start -1  --end 1
done


####### quality control

for i in mock Mo6hpi Mo12hpi Mo24hpi  ; do
Rscript $script/spatial_obj_qc.r \
--rds $workdir/obj_albst/${i}.rds  \
--nGene.min 50 \
--nspot.min 10 \
-o $workdir/obj_qc  -p ${i}
done


####### data scale

for i in  mock Mo6hpi Mo12hpi Mo24hpi ; do
Rscript $script/spatial_obj_scale.r \
-r $workdir/obj_qc/${i}.rds \
-o $workdir/obj_scale \
-p ${i}
done



####### gene rate heatmap plot
for j in  NLR  ; do
for i in mock Mo6hpi Mo12hpi Mo24hpi ; do
Rscript $script/spatial_albst_rate_heatmap.r \
-r $workdir/obj_scale/${i}_scale.rds  \
-l $workdir/list/${j}.txt \
-o $workdir/heatmap_rate/${j}  -H 0.1  -W 0.2 -p ${i}_${j}_rate \
--cluster  --heatmap 30  --fit 5 --gnum 7.5  --miny 1.5   --num 20  \
--start -1 --end 1 
done
done



################## ALBST cooperation with monocle2
######## monocle2 object construction

for i in  mock Mo6hpi Mo12hpi Mo24hpi  ; do
Rscript $script/spatial_monocle_construction.r \
 -i $workdir/obj_qc/${i}.rds -p ${i} \
 --heatmap.gene 200 --heatmap.clusters 3 \
 -o $workdir/monocle/${i}
done



###### monocle2 heatmap plot:

for j in ion; do
for i in mock Mo6hpi Mo12hpi Mo24hpi ; do
Rscript $script/spatial_monocle_heatmap_list.r \
 -i $workdir/monocle/${i}/${i}.spatial_trajectory.rds -p ${i}_${j} \
 --heatmap.clusters 1 \
 --list $workdir/list/${j}.txt  \
 -o $workdir/heatmap_monocle/${j}  -H 1 
done
done





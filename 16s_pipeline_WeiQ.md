## 16s pipeline for CALM2004 project
written by Wei Qing <br>
updated by Cancan Qi

## Prepare data
#### 1. 将总fastq文件中的序列拆分到各个样本，并去除barcode和primer
```bash
perl /data/qingwei/script/split_fq_according_barcode_parallel.pl -f <forward_seq> -r <reverse_seq> -m <metadata> -o <output_dir> -n <threads_num>
```
-f 正向序列文件 <br>
-r 反向序列文件 <br>
-m meta表，要包含#SampleID、Forward_Barcode、Forward_Primer、Reverse_Barcode、Reverse_Primer等信息 <br>
-o 结果目录 <br>
-n 线程数 <br>

#### 2. 将各样本用两套引物扩增的序列合并为一个fastq文件（在拆分到各样本的fastq文件夹中运行）
```bash
ls *.fastq > list.txt
sed -i 's/[A|B]_R1.fastq//g' list.txt
sed -i 's/[A|B]_R2.fastq//g' list.txt
uniq list.txt > uniq_list.txt
mkdir merged
for i in `tail -n+1 uniq_list.txt`
do
cat ${i}A_R1.fastq ${i}B_R1.fastq > merged/${i}_R1.fastq
cat ${i}A_R2.fastq ${i}B_R2.fastq > merged/${i}_R2.fastq
done
```
因测序未测通，后续将使用各样本所有的前向序列（*_R1.fastq）进行后续分析

#### 3. 构建样本目录清单文件
文件格式（制表符分割的文本文件，共两列）：<br>
sample-id	absolute-filepath <br>
sample-1	前向序列绝对路径：xxx/xxx/xxx/sample1_R1.fastq <br>
sample-2	xxx/xxx/xxx/sample2_R1.fastq <br>

在保存merged file的目录下
```bash
ls >filelist.txt
```
然后用R将该表做成上述格式
```R
setwd("./forward")
tb<-read.table("../filelist.txt")
tb$sampleid<-sub("_[^_]*$", "", tb$V1)
tb$sampleid[grep("NTC.",tb$sampleid)]<-paste0(tb$sampleid[grep("NTC.",tb$sampleid)],".")
tb$sampleid<-gsub('.{1}$', '', tb$sampleid)
tb$`absolute-filepath`<-paste0(getwd(),"/",tb$V1)
tb1<-tb[-dim(tb)[1],c("sampleid","absolute-filepath")] ## according to your own data
# tb1<-tb1[-grep("R2",tb1$`absolute-filepath`),] ## pleasse run this if the folder contains both R1 and R2
write.table(tb1,file="../file_dir_list.txt",quote=F,sep="\t",row.names=F)
```

## qiime2 mapping

#### 1.启动qiime2
```bash
source ~/miniconda3/bin/activate
conda activate qiime2-2022.2
```

#### 2.导入各样本前向fastq文件 并输出qzv表用于检查
```bash
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path 样本目录清单文件 \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2
  
qiime demux summarize \
  --i-data single-end-demux.qza \
  --o-visualization single-end-demux.qzv
```
随后即可以在 https://view.qiime2.org/ 上查看Demultiplexed sequence counts summary

#### 3.使用dada2去噪
```bash
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-demux.qza \
  --p-trunc-len 223 \
  --p-n-threads 线程数 \
  --o-table table.qza（此为ASV表） \
  --o-representative-sequences rep-seqs.qza（此为代表序列） \
  --o-denoising-stats denoising-stats.qza
  
## for data visulization
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file ## to be confirmed
```

#### 4. 基于代表性序列构建物种进化树
```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza \
  --p-n-threads 线程数
```

#### 5. 将代表性序列基于silva数据库进行物种分类
```bash
qiime feature-classifier classify-sklearn \
  --i-classifier /data/qingwei/CALM2004/classifier/silva/silva_full_length_classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza（此为物种分类结果文件） \
  --p-n-jobs 线程数
```

#### 6.删除所有分类为“Unassigned”、“d__Bacteria”和“d__Eukaryota”的序列
```bash
qiime taxa filter-table \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude "Unassigned" \
  --o-filtered-table filtered_table1.qza
  
qiime taxa filter-table \
  --i-table filtered_table1.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude "d__Bacteria;__;__;__;__;__;__" \
  --o-filtered-table filtered_table2.qza \
  --p-mode "contains"
  
qiime taxa filter-table \
  --i-table filtered_table2.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude "d__Eukaryota;__;__;__;__;__;__" \
  --o-filtered-table filtered_table3.qza \
  --p-mode "exact"
```

#### 7. 输出文件（biom.txt）可进行下游分析
```bash
qiime feature-table transpose \
  --i-table table.qza \
  --o-transposed-feature-table transposed-table.qza
  
qiime metadata tabulate \
  --m-input-file rep-seqs.qza \
  --m-input-file taxonomy.qza \
  --m-input-file transposed-table.qza \
  --o-visualization biom.qzv
  
qiime tools export \
  --input-path biom.qzv 
  --output-path biom


```



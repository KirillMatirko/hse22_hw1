# hse22_hw1

### 1. Создаём ссылки на файлы

```bash
ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq
```

### 2. Выбираем случайные чтения

```bash
seqtk sample -s927 oil_R1.fastq 5000000 > sub1.fq
seqtk sample -s927 oil_R2.fastq 5000000 > sub2.fq
seqtk sample -s927 oilMP_S4_L001_R1_001.fastq 1500000 > subMP1.fq
seqtk sample -s927 oilMP_S4_L001_R2_001.fastq 1500000 > subMP2.fq
```

### 3. Оценка качества чтений

```bash
mkdir fastqc
ls sub* subMP* | xargs -P 2 -tI{} fastqc -o fastqc {}
mkdir multiqc
multiqc -o multiqc fastqc
```
### 4. Скриншоты с multiqc
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/gen_stat.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/seq_counts.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/mean_scores.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/per_seq_scores.png)

### 5. Подрезка чтений и получение по ним статистики

```bash
platanus_trim sub
platanus_internal_trim subMP*
mkdir multiqc_trimmed
mkdir fastqc_trimmed
ls sub* matep*| xargs -tI{} fastqc -o fastqc_trimmed {}
multiqc -o multiqc_trimmed fastqc_trimmed
```

### 6. Удаление исходных чтений

```bash
rm sub1.fq sub2.fq subMP1.fq subMP2.fq
```

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

### 6. Скриншоты для подрезанных чтений

![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/gen_stat_trimmed.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/seq_counts_trimmed.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/mean_scores_trimmed.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/per_seq_scores_trimmed.png)

### 7. Удаление исходных чтений

```bash
rm sub1.fq sub2.fq subMP1.fq subMP2.fq
```

### 8. Сбор контигов

```bash
screen platanus assemble -o Poil -f sub1.fq.trimmed sub2.fq.trimmed 2
```

### 9. Анализ контигов

#### a) Питоновская функция, которую потом используем и для скаффолдов

```python
def analysis(file,label):
  lengths=[]
  l=0
  max_seq=''
  seq=''

  for line in file:
    if line[0]=='>':
      if l > 0:
        if len(seq) > len(max_seq):
          max_seq=seq
        seq=''
        lengths.append(l)
        l=0
      else:
        pass
    else:
       l+=len(line.strip())
       seq += line.strip()
  
  if len(seq) > len(max_seq):
          max_seq=seq
  
  lengths.append(l)
  lengths.sort(reverse=True)
  num_of_contings=len(lengths)
  total_length=sum(lengths)
  max_length=lengths[0]
  N50=0

  criterion=0
  for el in lengths:
    if criterion < total_length/2:
      criterion += el
    else:
      N50=el
      break

  print(f'{label}\n\
Число: {num_of_contings}\n\
Суммарная длина: {total_length}\n\
Максимальная длина: {max_length}\n\
N50: {N50}')
  return max_seq
```

#### б) Анализ данных для контигов

```python
output_contig = analysis(open('Poil_contig.fa', 'r'),'Контиги')
```
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/contig_analysis.png)

### 10. Сбор скаффолдов

```bash
screen platanus scaffold -o Poil -c Poil_contig.fa -IP1 sub1.fq.trimmed sub2.fq.trimmed -OP2 subMP1.fq.trimmed subMP2.fq.trimmed 2
```
### 11. Анализ данных для скаффолдов

```python
output_scaffold = analysis(open('Poil_scaffold.fa', 'r'),'Скаффолды')
```
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/scaffold_analysis.png)

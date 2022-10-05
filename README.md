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
seqtk sample -s927 oil_R1.fastq 5000000 > sample1.fastq
seqtk sample -s927 oil_R2.fastq 5000000 > sample2.fastq
seqtk sample -s927 oilMP_S4_L001_R1_001.fastq 1500000 > matched_sample1.fastq
seqtk sample -s927 oilMP_S4_L001_R2_001.fastq 1500000 > matched_sample2.fastq
```

### 3. Оценка качества чтений

```bash
mkdir fastqc
mkdir multiqc
fastqc -o fastqc sample1.fastq sample2.fastq matched_sample1.fastq matched_sample2.fastq
multiqc -o multiqc fastqc
```
### 4. Скриншоты с multiqc
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/gen_stat.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/seq_counts.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/mean_scores.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/per_seq_scores.png)

### 5. Подрезка чтений и получение по ним статистики

```bash
platanus_trim sample1.fastq sample2.fastq
platanus_internal_trim matched_sample1.fastq matched_sample2.fastq
mkdir multiqc_trimmed
mkdir fastqc_trimmed
fastqc -o fastqc_trimmed sample1.fastq.trimmed sample2.fastq.trimmed matched_sample1.fastq.int_trimmed matched_sample2.fastq.int_trimmed
multiqc -o multiqc_trimmed fastqc_trimmed
```

### 6. Скриншоты для подрезанных чтений

![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/gen_stat_trimmed.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/seq_counts_trimmed.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/mean_scores_trimmed.png)
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/per_seq_scores_trimmed.png)

### 7. Удаление исходных чтений

```bash
rm sample1.fastq sample2.fastq matched_sample1.fastq matched_sample2.fastq
```

### 8. Сбор контигов

```bash
screen platanus assemble -o Poil -f sub1.fq.trimmed sub2.fq.trimmed 2
```

### 9. Анализ контигов

Ссылка на Colab https://colab.research.google.com/drive/1K6Rqg3CRt1OJiBYv9Wfe1fKg7iyU3qST?usp=sharing

#### a) Питоновская функция, которую потом используем и для скаффолдов

```python
import re

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

### 12. Считаем число гэпов для самого длинного скаффолда

```python
length_of_gaps = output_scaffold.count('N')
output_scaffold_new = re.sub(r'N{2,}', 'N', output_scaffold)
num_of_gaps = output_scaffold_new.count('N')
print(f'Длина гэпов: {length_of_gaps}\n\
Число гэпов: {num_of_gaps}')
```
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/gaps.png)

### 13. Уменьшаем число гэпов

```bash
screen platanus gap_close -o Poil -c Poil_scaffold.fa -IP1 sub1.fq.trimmed sub2.fq.trimmed -OP2 subMP1.fq.trimmed subMP2.fq.trimmed 2
```

### 14. Считаем число гэпов для самого длинного скаффолда в укороченной версии

```python
output_gapClosed = analysis(open('Poil_gapClosed.fa', 'r'),'Скаффолды с укороченными гэпами')
```
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/scaffold_truncated_analysis.png)

```python
length_of_gaps2 = output_gapClosed.count('N')
output_gapClosed_new = re.sub(r'N{2,}', 'N', output_gapClosed)
num_of_gaps2 = output_gapClosed_new.count('N')
print(f'Укороченная версия\n\
Длина гэпов: {length_of_gaps2}\n\
Число гэпов: {num_of_gaps2}')
```
![](https://github.com/KirillMatirko/hse22_hw1/blob/main/pics/gaps_truncated.png)

### 15. Записываем самые длинные скаффолды

```python
with open('longest_scaffold.fa', 'w') as f:
  f.write(output_scaffold)
with open('longest_truncated_scaffold.fa', 'w') as f:
  f.write(output_gapClosed)
```

### 16. Удаляем подрезанные чтения

```bash
rm sub*
```

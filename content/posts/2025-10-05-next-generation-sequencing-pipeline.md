---
title: "构建完整的NGS数据分析流程"
publishedAt: "2025-10-05"
updatedAt: "2025-10-05"
excerpt: "详细介绍从原始测序数据到变异检测的完整NGS分析流程，包括质量控制、比对、变异 calling等步骤"
tags: ["NGS", "测序", "生物信息学", "流程", "Python"]
authors: ["陈七", "周八"]
coverImage: "https://images.unsplash.com/photo-1581093450021-4a7360e9a6b5?w=1200&h=630&fit=crop&crop=center&auto=format&q=80"
---

# 构建完整的NGS数据分析流程

随着二代测序（Next-Generation Sequencing, NGS）技术的普及，构建高效、准确的数据分析流程变得至关重要。本文将详细介绍从原始测序数据到最终变异检测的完整分析流程。

## 1. 数据质量控制

### 1.1 FastQC质量检查

```bash
# 运行FastQC检查原始数据质量
fastqc sample_R1.fastq.gz sample_R2.fastq.gz -o fastqc_output/

# 使用MultiQC汇总多个样本的质控结果
multiqc fastqc_output/ -o multiqc_output/
```

### 1.2 Trimmomatic数据清理

```bash
java -jar trimmomatic-0.39.jar PE \
  sample_R1.fastq.gz sample_R2.fastq.gz \
  sample_R1_trimmed.fastq.gz sample_R1_unpaired.fastq.gz \
  sample_R2_trimmed.fastq.gz sample_R2_unpaired.fastq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

## 2. 序列比对

### 2.1 BWA比对

```bash
# 构建参考基因组索引
bwa index reference.fa

# 执行比对
bwa mem -t 8 reference.fa \
  sample_R1_trimmed.fastq.gz sample_R2_trimmed.fastq.gz \
  > sample.sam

# 转换为BAM格式并排序
samtools view -bS sample.sam | samtools sort -o sample_sorted.bam
samtools index sample_sorted.bam
```

### 2.2 Python自动化比对流程

```python
import subprocess
import os
from pathlib import Path

class NGSAnalysisPipeline:
    def __init__(self, config):
        self.config = config
        self.sample_dir = Path(config['sample_dir'])
        self.output_dir = Path(config['output_dir'])
        self.reference = Path(config['reference'])

    def run_bwa_alignment(self, sample_name):
        """执行BWA比对"""
        r1 = self.sample_dir / f"{sample_name}_R1_trimmed.fastq.gz"
        r2 = self.sample_dir / f"{sample_name}_R2_trimmed.fastq.gz"

        # 创建输出目录
        align_dir = self.output_dir / "alignment" / sample_name
        align_dir.mkdir(parents=True, exist_ok=True)

        # BWA比对
        sam_file = align_dir / f"{sample_name}.sam"
        cmd = f"bwa mem -t 8 {self.reference} {r1} {r2} > {sam_file}"
        subprocess.run(cmd, shell=True, check=True)

        # 转换为BAM并排序
        bam_file = align_dir / f"{sample_name}_sorted.bam"
        subprocess.run(f"samtools view -bS {sam_file} | samtools sort -o {bam_file}",
                      shell=True, check=True)
        subprocess.run(f"samtools index {bam_file}", shell=True, check=True)

        return bam_file
```

## 3. 比对后处理

### 3.1 Picard工具处理

```bash
# 添加Read组信息
java -jar picard.jar AddOrReplaceReadGroups \
  I=sample_sorted.bam \
  O=sample_rg.bam \
  RGID=sample1 \
  RGLB=lib1 \
  RGPL=ILLUMINA \
  RGPU=unit1 \
  RGSM=sample1

# 标记重复序列
java -jar picard.jar MarkDuplicates \
  I=sample_rg.bam \
  O=sample_dedup.bam \
  M=metrics.txt
```

### 3.2 GATK最佳实践流程

```python
def run_gatk_base_recalibration(self, bam_file, sample_name):
    """运行GATK碱基质量重新校准"""
    gatk_dir = self.output_dir / "gatk" / sample_name
    gatk_dir.mkdir(parents=True, exist_ok=True)

    # 已知变异位点
    known_sites = [
        self.config['dbsnp'],
        self.config['mills'],
        self.config['indels']
    ]

    # BaseRecalibrator
    recal_table = gatk_dir / f"{sample_name}_recal.table"
    cmd = (f"gatk BaseRecalibrator "
           f"-I {bam_file} "
           f"-R {self.reference} "
           f"{' '.join([f'--known-sites {site}' for site in known_sites])} "
           f"-O {recal_table}")
    subprocess.run(cmd, shell=True, check=True)

    # ApplyBQSR
    recal_bam = gatk_dir / f"{sample_name}_recal.bam"
    cmd = (f"gatk ApplyBQSR "
           f"-I {bam_file} "
           f"-R {self.reference} "
           f"--bqsr-recal-file {recal_table} "
           f"-O {recal_bam}")
    subprocess.run(cmd, shell=True, check=True)

    return recal_bam
```

## 4. 变异检测

### 4.1 GATK HaplotypeCaller

```bash
# 单样本变异检测
gatk HaplotypeCaller \
  -R reference.fa \
  -I sample_recal.bam \
  -O sample_raw.vcf.gz

# 联合基因分型
gatk GenomicsDBImport \
  --sample-name-map sample_map.txt \
  --genomicsdb-workspace-path my_database \
  --reference-fasta reference.fa

gatk GenotypeGVCFs \
  -R reference.fa \
  -V gendb://my_database \
  -O cohort.vcf.gz
```

### 4.2 Python变异过滤

```python
import pysam

def filter_variants(vcf_file, output_file, min_quality=30, min_depth=10):
    """过滤低质量变异"""
    vcf_reader = pysam.VariantFile(vcf_file)
    vcf_writer = pysam.VariantFile(output_file, 'w', header=vcf_reader.header)

    for record in vcf_reader:
        # 过滤条件
        if record.qual < min_quality:
            continue

        # 检查测序深度
        if any(sample.get('DP', 0) < min_depth for sample in record.samples):
            continue

        # 写入过滤后的变异
        vcf_writer.write(record)

    vcf_reader.close()
    vcf_writer.close()

# 变异注释
def annotate_variants(vcf_file, annotation_db):
    """使用VEP注释变异"""
    cmd = f"vep -i {vcf_file} -o {vcf_file}.annotated.vcf --cache --offline"
    subprocess.run(cmd, shell=True, check=True)
```

## 5. 变异功能分析

### 5.1 SnpEff功能预测

```bash
# 使用SnpEff进行功能注释
java -jar snpEff.jar -v GRCh38.99 cohort_filtered.vcf > cohort_annotated.vcf

# 生成注释统计报告
java -jar snpEff.jar -stats snpEff_stats.html GRCh38.99 cohort_filtered.vcf
```

### 5.2 Python变异分析

```python
import pandas as pd

class VariantAnalyzer:
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.variants_df = self._parse_vcf()

    def _parse_vcf(self):
        """解析VCF文件"""
        variants = []
        vcf_reader = pysam.VariantFile(self.vcf_file)

        for record in vcf_reader:
            variant_info = {
                'chrom': record.chrom,
                'pos': record.pos,
                'ref': record.ref,
                'alt': record.alts[0] if record.alts else '.',
                'qual': record.qual,
                'filter': record.filter,
                'type': self._get_variant_type(record),
                'impact': record.info.get('ANN', [''])[0].split('|')[2] if 'ANN' in record.info else 'UNKNOWN'
            }
            variants.append(variant_info)

        return pd.DataFrame(variants)

    def _get_variant_type(self, record):
        """确定变异类型"""
        if len(record.ref) == 1 and len(record.alts) == 1 and len(record.alts[0]) == 1:
            return 'SNP'
        elif len(record.ref) < len(record.alts[0]):
            return 'INS'
        elif len(record.ref) > len(record.alts[0]):
            return 'DEL'
        else:
            return 'COMPLEX'

    def get_variant_statistics(self):
        """获取变异统计信息"""
        stats = {
            'total_variants': len(self.variants_df),
            'snp_count': len(self.variants_df[self.variants_df['type'] == 'SNP']),
            'indel_count': len(self.variants_df[self.variants_df['type'].isin(['INS', 'DEL'])]),
            'high_impact': len(self.variants_df[self.variants_df['impact'] == 'HIGH']),
            'moderate_impact': len(self.variants_df[self.variants_df['impact'] == 'MODERATE'])
        }
        return stats
```

## 6. 质量评估与报告

### 6.1 使用Qualimap评估比对质量

```bash
# 运行Qualimap
qualimap bamqc -bam sample_final.bam -outdir qualimap_report

# 生成HTML报告
qualimap multi-bamqc -d bam_list.txt -outdir multi_qualimap_report
```

### 6.2 Python自动化报告生成

```python
import json
from datetime import datetime

class NGSReportGenerator:
    def __init__(self, analysis_results):
        self.results = analysis_results
        self.report_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    def generate_html_report(self, output_file):
        """生成HTML格式报告"""
        html_template = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>NGS分析报告</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background-color: #f0f0f0; padding: 20px; }}
                .section {{ margin: 20px 0; }}
                .metric {{ display: inline-block; margin: 10px; padding: 10px; border: 1px solid #ddd; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>NGS数据分析报告</h1>
                <p>生成时间: {self.report_date}</p>
            </div>

            <div class="section">
                <h2>样本信息</h2>
                <p>样本数量: {self.results['sample_count']}</p>
                <p>测序平台: {self.results['sequencing_platform']}</p>
            </div>

            <div class="section">
                <h2>比对统计</h2>
                <div class="metric">
                    <h3>平均比对率</h3>
                    <p>{self.results['alignment_rate']:.2f}%</p>
                </div>
                <div class="metric">
                    <h3>平均覆盖深度</h3>
                    <p>{self.results['mean_coverage']:.1f}x</p>
                </div>
            </div>

            <div class="section">
                <h2>变异统计</h2>
                <table>
                    <tr><th>变异类型</th><th>数量</th></tr>
                    <tr><td>SNP</td><td>{self.results['variant_stats']['snp_count']}</td></tr>
                    <tr><td>插入/缺失</td><td>{self.results['variant_stats']['indel_count']}</td></tr>
                    <tr><td>高影响变异</td><td>{self.results['variant_stats']['high_impact']}</td></tr>
                </table>
            </div>
        </body>
        </html>
        """

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_template)
```

## 7. 完整流程自动化

### 7.1 使用Snakemake构建流程

```python
# Snakefile
configfile: "config.yaml"

rule all:
    input:
        "reports/final_report.html"

rule fastqc:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        html = "qc/fastqc/{sample}_R1_fastqc.html"
    shell:
        "fastqc {input.r1} {input.r2} -o qc/fastqc/"

rule trimmomatic:
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
    output:
        r1 = "trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "trimmed/{sample}_R2_trimmed.fastq.gz"
    shell:
        "trimmomatic PE {input.r1} {input.r2} {output.r1} /dev/null {output.r2} /dev/null "
        "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule bwa_align:
    input:
        r1 = "trimmed/{sample}_R1_trimmed.fastq.gz",
        r2 = "trimmed/{sample}_R2_trimmed.fastq.gz",
        ref = config['reference']
    output:
        bam = "aligned/{sample}_sorted.bam"
    shell:
        "bwa mem -t 8 {input.ref} {input.r1} {input.r2} | "
        "samtools sort -o {output.bam} && "
        "samtools index {output.bam}"

rule generate_report:
    input:
        bams = expand("aligned/{sample}_sorted.bam", sample=config['samples']),
        vcf = "variants/cohort_filtered.vcf"
    output:
        "reports/final_report.html"
    script:
        "scripts/generate_report.py"
```

## 8. 性能优化与最佳实践

### 8.1 并行化处理

```python
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

def parallel_bwa_alignment(samples, max_workers=None):
    """并行化BWA比对"""
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = []
        for sample in samples:
            future = executor.submit(run_bwa_alignment, sample)
            futures.append(future)

        results = []
        for future in futures:
            results.append(future.result())

    return results
```

### 8.2 内存优化

```python
def memory_efficient_vcf_processing(vcf_file, chunk_size=10000):
    """内存高效的VCF文件处理"""
    vcf_reader = pysam.VariantFile(vcf_file)

    # 分批处理
    batch = []
    for record in vcf_reader:
        batch.append(record)

        if len(batch) >= chunk_size:
            process_variant_batch(batch)
            batch = []

    # 处理最后一批
    if batch:
        process_variant_batch(batch)

    vcf_reader.close()
```

## 总结

本文详细介绍了完整的NGS数据分析流程，从原始数据质量控制到最终变异检测和功能注释。通过Python自动化和流程管理工具，可以构建高效、可重复的分析流程。随着测序技术的不断发展，优化和标准化分析流程将变得越来越重要。

## 参考资料

1. GATK Best Practices: https://gatk.broadinstitute.org/hc/en-us/articles/360035531192
2. SAMtools Documentation: https://www.htslib.org/doc/samtools.html
3. Picard Tools: https://broadinstitute.github.io/picard/
4. Snakemake Workflow Management: https://snakemake.readthedocs.io/
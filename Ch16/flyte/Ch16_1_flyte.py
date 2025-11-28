#!/usr/bin/env python3
"""
Bioinformatics Workflow with Flyte - Command Line Version

A complete genomics pipeline for quality control, alignment, and variant calling.
Designed to be run from the command line using pyflyte.

Usage:
    # First, generate sample files 
    python Ch16_1_flyte.py
    # Run locally
    pyflyte run Ch16_1_flyte.py genomics_pipeline --fastq_file sample.fastq --reference_genome reference.fa --sample_name "my_sample"
    
    # Run with Flyte sandbox
    flytectl demo start
    pyflyte --config ~/.flyte/config-sandbox.yaml register Ch16_1_flyte.py
    pyflyte --config ~/.flyte/config-sandbox.yaml run --remote Ch16_1_flyte.py genomics_pipeline --fastq_file sample.fastq --reference_genome reference.fa
    
Requirements:
    pip install "flytekit[all]" 

Author: Flyte Genomics Pipeline
Version: 1.0
"""

# Import Libraries
import os
import tempfile
import random
import argparse
import warnings
from typing import Tuple
from dataclasses import dataclass

# Suppress specific warnings from flytekit/click
warnings.filterwarnings("ignore", message=".*parameter.*is used more than once.*", category=UserWarning)
warnings.filterwarnings("ignore", message=".*parameter -i is used more than once.*", category=UserWarning)

from flytekit import task, workflow, Resources
from flytekit.types.file import FlyteFile
from flytekit.types.directory import FlyteDirectory

# Dataclass Definitions
@dataclass
class QCMetrics:
    """Data class to hold quality control metrics"""
    total_reads: int
    mean_quality: float
    gc_content: float
    duplication_rate: float

    def __str__(self):
        return (f"QCMetrics(reads={self.total_reads}, quality={self.mean_quality:.2f}, "
                f"gc={self.gc_content:.2f}%, dup_rate={self.duplication_rate:.2f}%)")


@dataclass
class AlignmentMetrics:
    """Data class to hold alignment metrics"""
    total_reads: int
    mapped_reads: int
    mapping_rate: float
    mean_coverage: float

    def __str__(self):
        return (f"AlignmentMetrics(total={self.total_reads}, mapped={self.mapped_reads}, "
                f"rate={self.mapping_rate:.2f}%, coverage={self.mean_coverage:.2f}x)")

# Task Definitions
@task(
    requests=Resources(cpu="1", mem="1Gi"),
    cache_version="1.0"
)
def quality_control(fastq_file: FlyteFile) -> Tuple[FlyteDirectory, QCMetrics]:
    """
    Perform quality control analysis on FASTQ files
    
    Args:
        fastq_file: Input FASTQ file for analysis
        
    Returns:
        Tuple of QC output directory and metrics
    """
    print("🔬 Starting Quality Control Analysis...")
    
    # Create output directory
    output_dir = tempfile.mkdtemp(prefix="qc_output_")
    
    # Download the input file to analyze
    fastq_path = fastq_file.download()
    print(f"   Processing file: {fastq_path}")
    
    # Parse FASTQ file to extract real metrics
    total_reads = 0
    total_quality = 0
    total_bases = 0
    gc_count = 0
    
    try:
        with open(fastq_path, 'r') as f:
            lines = f.readlines()
            
        # Process FASTQ format (4 lines per read)
        for i in range(0, len(lines), 4):
            if i + 3 < len(lines) and lines[i].startswith('@'):
                total_reads += 1
                
                # Sequence line
                sequence = lines[i + 1].strip()
                total_bases += len(sequence)
                gc_count += sequence.count('G') + sequence.count('C')
                
                # Quality line
                quality_line = lines[i + 3].strip()
                # Convert ASCII quality scores to numeric (Phred+33)
                quality_scores = [ord(c) - 33 for c in quality_line]
                total_quality += sum(quality_scores)
        
        # Calculate metrics
        mean_quality = total_quality / total_bases if total_bases > 0 else 30
        gc_content = (gc_count / total_bases * 100) if total_bases > 0 else 45
        
    except Exception as e:
        print(f"   Warning: Error parsing FASTQ file: {e}")
        # Use default values if parsing fails
        total_reads = 100000
        mean_quality = 32.5
        gc_content = 42.1
    
    qc_metrics = QCMetrics(
        total_reads=total_reads,
        mean_quality=mean_quality,
        gc_content=gc_content,
        duplication_rate=random.uniform(10, 20)  # Simulated
    )
    
    # Create detailed QC report files
    report_file = os.path.join(output_dir, "fastqc_report.html")
    with open(report_file, 'w') as f:
        f.write(f"""
<!DOCTYPE html>
<html>
<head>
    <title>FastQC Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        table {{ border-collapse: collapse; width: 100%; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .pass {{ color: green; }}
        .warn {{ color: orange; }}
        .fail {{ color: red; }}
    </style>
</head>
<body>
    <h1>FastQC Quality Control Report</h1>
    <h2>Basic Statistics</h2>
    <table>
        <tr><th>Metric</th><th>Value</th><th>Status</th></tr>
        <tr><td>Total Reads</td><td>{qc_metrics.total_reads:,}</td>
            <td class="{'pass' if qc_metrics.total_reads > 1000 else 'warn'}">
                {'PASS' if qc_metrics.total_reads > 1000 else 'WARN'}
            </td></tr>
        <tr><td>Mean Quality Score</td><td>{qc_metrics.mean_quality:.2f}</td>
            <td class="{'pass' if qc_metrics.mean_quality > 25 else 'warn' if qc_metrics.mean_quality > 15 else 'fail'}">
                {'PASS' if qc_metrics.mean_quality > 25 else 'WARN' if qc_metrics.mean_quality > 15 else 'FAIL'}
            </td></tr>
        <tr><td>GC Content</td><td>{qc_metrics.gc_content:.2f}%</td>
            <td class="{'pass' if 30 <= qc_metrics.gc_content <= 70 else 'warn'}">
                {'PASS' if 30 <= qc_metrics.gc_content <= 70 else 'WARN'}
            </td></tr>
        <tr><td>Duplication Rate</td><td>{qc_metrics.duplication_rate:.2f}%</td>
            <td class="{'pass' if qc_metrics.duplication_rate < 30 else 'warn'}">
                {'PASS' if qc_metrics.duplication_rate < 30 else 'WARN'}
            </td></tr>
    </table>
    
    <h2>Quality Assessment</h2>
    <p><strong>Overall Status:</strong> 
       <span class="{'pass' if qc_metrics.mean_quality > 25 and 30 <= qc_metrics.gc_content <= 70 else 'warn'}">
           {'PASS - Good quality data' if qc_metrics.mean_quality > 25 and 30 <= qc_metrics.gc_content <= 70 else 'REVIEW - Check quality metrics'}
       </span>
    </p>
    
    <h2>Recommendations</h2>
    <ul>
        {'<li>Data quality is acceptable for downstream analysis</li>' if qc_metrics.mean_quality > 25 else '<li>Consider quality filtering or trimming</li>'}
        {'<li>GC content is within normal range</li>' if 30 <= qc_metrics.gc_content <= 70 else '<li>Unusual GC content - check for contamination</li>'}
        {'<li>Duplication levels are acceptable</li>' if qc_metrics.duplication_rate < 30 else '<li>High duplication detected - consider deduplication</li>'}
    </ul>
</body>
</html>
        """)
    
    # Create summary stats file
    stats_file = os.path.join(output_dir, "qc_summary.txt")
    with open(stats_file, 'w') as f:
        f.write(f"Quality Control Summary\n")
        f.write(f"======================\n")
        f.write(f"Total Reads: {qc_metrics.total_reads:,}\n")
        f.write(f"Mean Quality: {qc_metrics.mean_quality:.2f}\n")
        f.write(f"GC Content: {qc_metrics.gc_content:.2f}%\n")
        f.write(f"Duplication Rate: {qc_metrics.duplication_rate:.2f}%\n")
    
    print(f"   ✅ QC Complete: {qc_metrics}")
    print(f"   📁 Output directory: {output_dir}")
    
    return FlyteDirectory(path=output_dir), qc_metrics


@task(
    requests=Resources(cpu="2", mem="2Gi"),
    cache_version="1.0"
)
def align_reads(
    fastq_file: FlyteFile, 
    reference_genome: FlyteFile,
    sample_name: str = "sample"
) -> Tuple[FlyteFile, AlignmentMetrics]:
    """
    Align reads to reference genome (simulated for local execution)
    
    Args:
        fastq_file: Input FASTQ file
        reference_genome: Reference genome FASTA file
        sample_name: Sample identifier
        
    Returns:
        Tuple of aligned BAM file and alignment metrics
    """
    print("🧬 Starting Read Alignment...")
    
    # Download input files
    fastq_path = fastq_file.download()
    ref_path = reference_genome.download()
    print(f"   FASTQ: {fastq_path}")
    print(f"   Reference: {ref_path}")
    
    # Count reads from FASTQ file
    try:
        with open(fastq_path, 'r') as f:
            lines = f.readlines()
            total_reads = len([line for line in lines if line.startswith('@')])
    except Exception:
        total_reads = 100000  # Default for simulation
    
    # Create a simulated BAM file (text format for simplicity)
    sorted_bam = f"{sample_name}_sorted.bam"
    with open(sorted_bam, 'w') as f:
        f.write(f"# Simulated BAM file for {sample_name}\n")
        f.write(f"# Total reads processed: {total_reads}\n")
        f.write(f"# Generated by Flyte Genomics Pipeline\n")
        f.write("# SAM Header\n")
        f.write("@HD\tVN:1.0\tSO:coordinate\n")
        f.write("@SQ\tSN:chr1\tLN:248956422\n")
        f.write("@SQ\tSN:chr2\tLN:242193529\n")
        f.write("@RG\tID:sample\tSM:" + sample_name + "\n")
        f.write("# Alignment records (first 20 shown)\n")
        
        # Add some mock alignment records
        chromosomes = ["chr1", "chr2"]
        bases = ["A", "T", "G", "C"]
        for i in range(min(20, total_reads)):
            chrom = random.choice(chromosomes)
            pos = random.randint(1000, 1000000)
            mapq = random.randint(20, 60)
            seq = ''.join(random.choice(bases) for _ in range(36))
            qual = 'I' * 36  # High quality
            f.write(f"read_{i+1}\t0\t{chrom}\t{pos}\t{mapq}\t36M\t*\t0\t0\t{seq}\t{qual}\n")
    
    # Calculate simulated alignment metrics
    mapping_rate = random.uniform(85, 95)  # 85-95% mapping rate
    mapped_reads = int(total_reads * (mapping_rate / 100))
    
    alignment_metrics = AlignmentMetrics(
        total_reads=total_reads,
        mapped_reads=mapped_reads,
        mapping_rate=mapping_rate,
        mean_coverage=random.uniform(25, 35)
    )
    
    print(f"   ✅ Alignment Complete: {alignment_metrics}")
    print(f"   📄 BAM file: {sorted_bam}")
    
    return FlyteFile(path=sorted_bam), alignment_metrics


@task(
    requests=Resources(cpu="1", mem="1Gi"),
    cache_version="1.0"
)
def call_variants(
    bam_file: FlyteFile, 
    reference_genome: FlyteFile,
    sample_name: str = "sample"
) -> FlyteFile:
    """
    Call variants from aligned reads (simulated)
    
    Args:
        bam_file: Aligned BAM file
        reference_genome: Reference genome FASTA file
        sample_name: Sample identifier
        
    Returns:
        VCF file with called variants
    """
    print("🔍 Starting Variant Calling...")
    
    # Download input files
    bam_path = bam_file.download()
    ref_path = reference_genome.download()
    print(f"   BAM: {bam_path}")
    print(f"   Reference: {ref_path}")
    
    vcf_file = f"{sample_name}_variants.vcf"
    
    # Create a realistic VCF file with mock variants
    with open(vcf_file, 'w') as f:
        # VCF header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=FlyteGenomicsPipeline\n")
        f.write("##reference=reference.fa\n")
        f.write(f"##sampleName={sample_name}\n")
        f.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n")
        f.write("##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n")
        f.write("##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant Type\">\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample_name + "\n")
        
        # Generate realistic variants
        chromosomes = ["chr1", "chr2"]
        bases = ["A", "T", "G", "C"]
        variant_types = ["SNP"] * 8 + ["INDEL"] * 2  # 80% SNPs, 20% INDELs
        
        num_variants = random.randint(100, 300)
        for i in range(num_variants):
            chrom = random.choice(chromosomes)
            pos = random.randint(10000, 900000)
            variant_type = random.choice(variant_types)
            
            if variant_type == "SNP":
                ref = random.choice(bases)
                alt = random.choice([b for b in bases if b != ref])
                info_type = "SNP"
            else:  # INDEL
                if random.random() < 0.6:  # 60% insertions
                    ref = random.choice(bases)
                    alt = ref + ''.join(random.choice(bases) for _ in range(random.randint(1, 3)))
                    info_type = "INS"
                else:  # 40% deletions
                    ref = ''.join(random.choice(bases) for _ in range(random.randint(2, 4)))
                    alt = ref[0]
                    info_type = "DEL"
            
            qual = random.randint(30, 90)
            depth = random.randint(15, 60)
            allele_freq = random.uniform(0.25, 0.75)
            alt_depth = int(depth * allele_freq)
            ref_depth = depth - alt_depth
            
            # Genotype based on allele frequency
            if allele_freq > 0.8:
                genotype = "1/1"  # Homozygous alternate
            elif allele_freq < 0.2:
                genotype = "0/0"  # Homozygous reference (shouldn't happen in variant calls)
            else:
                genotype = "0/1"  # Heterozygous
            
            filter_status = "PASS" if qual > 30 and depth > 10 else "LowQual"
            
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter_status}\t")
            f.write(f"DP={depth};AF={allele_freq:.3f};TYPE={info_type}\tGT:DP:AD\t")
            f.write(f"{genotype}:{depth}:{ref_depth},{alt_depth}\n")
    
    print(f"   ✅ Variant Calling Complete: {num_variants} variants found")
    print(f"   📄 VCF file: {vcf_file}")
    
    return FlyteFile(path=vcf_file)


@task(
    requests=Resources(cpu="1", mem="1Gi")
)
def generate_report(
    qc_metrics: QCMetrics,
    alignment_metrics: AlignmentMetrics,
    vcf_file: FlyteFile,
    sample_name: str = "sample"
) -> FlyteFile:
    """
    Generate comprehensive analysis report with visualizations
    
    Args:
        qc_metrics: Quality control metrics
        alignment_metrics: Alignment statistics  
        vcf_file: Variant call file
        sample_name: Sample identifier
        
    Returns:
        Comprehensive analysis report file
    """
    print("📊 Generating Analysis Report with Visualizations...")
    
    # Import standard library modules only
    import datetime
    import csv
    
    # Try to use matplotlib for visualizations (optional)
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend
        plt.style.use('default')
        visualizations_available = True
        print("   ✅ Matplotlib available - generating plots")
    except ImportError:
        # Try to install matplotlib
        try:
            import subprocess
            import sys
            print("   📦 Installing matplotlib for cluster execution...")
            subprocess.check_call([
                sys.executable, "-m", "pip", "install", 
                "--user", "--quiet", "matplotlib"
            ])
            import matplotlib.pyplot as plt
            import matplotlib
            matplotlib.use('Agg')
            plt.style.use('default')
            visualizations_available = True
            print("   ✅ Matplotlib installed and available")
        except Exception:
            visualizations_available = False
            print("   ⚠️ Matplotlib not available - skipping plots")
    
    # Download VCF file to analyze variants
    vcf_path = vcf_file.download()
    
    # Parse VCF file for variant statistics
    snp_count = 0
    indel_count = 0
    total_variants = 0
    high_quality_variants = 0
    
    with open(vcf_path, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                total_variants += 1
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    ref, alt, qual, filter_status, info = fields[3], fields[4], fields[5], fields[6], fields[7]
                    
                    # Count variant types
                    if len(ref) == 1 and len(alt) == 1:
                        snp_count += 1
                    else:
                        indel_count += 1
                    
                    # Count high quality variants
                    if filter_status == "PASS" and float(qual) > 30:
                        high_quality_variants += 1
    
    # Generate visualizations if matplotlib is available
    plot_files = []
    if visualizations_available:
        # Create 4-panel visualization
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # 1. Quality Metrics Bar Chart
        metrics = ['Quality Score', 'GC Content', 'Mapping Rate']
        values = [qc_metrics.mean_quality, qc_metrics.gc_content, alignment_metrics.mapping_rate]
        colors = ['skyblue', 'lightgreen', 'coral']
        
        bars = ax1.bar(metrics, values, color=colors)
        ax1.set_title('Key Quality Metrics', fontsize=14, fontweight='bold')
        ax1.set_ylabel('Score/Percentage')
        
        # Add value labels on bars
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                    f'{value:.1f}', ha='center', va='bottom')
        
        # 2. Read Mapping Pie Chart
        mapped = alignment_metrics.mapped_reads
        unmapped = alignment_metrics.total_reads - mapped
        
        ax2.pie([mapped, unmapped], 
                labels=[f'Mapped\n({mapped:,})', f'Unmapped\n({unmapped:,})'],
                autopct='%1.1f%%', 
                colors=['lightblue', 'lightcoral'],
                startangle=90)
        ax2.set_title('Read Mapping Distribution', fontsize=14, fontweight='bold')
        
        # 3. Quality Assessment Dashboard
        quality_categories = ['Read Quality', 'GC Content', 'Mapping', 'Coverage']
        quality_scores = [
            min(qc_metrics.mean_quality / 40 * 100, 100),  # Normalize to 100
            100 - abs(qc_metrics.gc_content - 50),  # Penalty for deviation from 50%
            alignment_metrics.mapping_rate,
            min(alignment_metrics.mean_coverage / 30 * 100, 100)  # Normalize to 100
        ]
        
        y_pos = range(len(quality_categories))
        bars = ax3.barh(y_pos, quality_scores, color=['green' if score > 80 else 'orange' if score > 60 else 'red' for score in quality_scores])
        ax3.set_yticks(y_pos)
        ax3.set_yticklabels(quality_categories)
        ax3.set_xlabel('Quality Score (%)')
        ax3.set_title('Quality Assessment Dashboard', fontsize=14, fontweight='bold')
        ax3.set_xlim(0, 100)
        
        # Add score labels
        for i, (bar, score) in enumerate(zip(bars, quality_scores)):
            width = bar.get_width()
            ax3.text(width + 1, bar.get_y() + bar.get_height()/2,
                    f'{score:.1f}%', ha='left', va='center')
        
        # 4. Variant Type Distribution
        if total_variants > 0:
            variant_labels = ['SNPs', 'INDELs']
            variant_counts = [snp_count, indel_count]
            variant_colors = ['lightsteelblue', 'lightsalmon']
            
            wedges, texts, autotexts = ax4.pie(variant_counts, labels=variant_labels, 
                                               autopct='%1.1f%%', colors=variant_colors,
                                               startangle=90)
            ax4.set_title(f'Variant Types\n(Total: {total_variants:,})', fontsize=14, fontweight='bold')
        else:
            ax4.text(0.5, 0.5, 'No Variants\nDetected', ha='center', va='center', 
                    transform=ax4.transAxes, fontsize=16)
            ax4.set_title('Variant Distribution', fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        # Save the main plot
        main_plot_file = f'{sample_name}_analysis_plots.png'
        plt.savefig(main_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        plot_files.append(main_plot_file)
        
        # Create additional individual plots
        
        # Quality trend plot
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        categories = ['Total Reads\n(thousands)', 'Mean Quality', 'GC Content (%)', 'Mapping Rate (%)']
        values = [qc_metrics.total_reads/1000, qc_metrics.mean_quality, qc_metrics.gc_content, alignment_metrics.mapping_rate]
        
        bars = ax.bar(categories, values, color=['steelblue', 'forestgreen', 'gold', 'crimson'])
        ax.set_title(f'Quality Metrics Summary - {sample_name}', fontsize=16, fontweight='bold')
        ax.set_ylabel('Value')
        
        # Add value labels
        for bar, value in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + max(values)*0.01,
                   f'{value:.1f}', ha='center', va='bottom', fontweight='bold')
        
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        
        quality_plot_file = f'{sample_name}_quality_summary.png'
        plt.savefig(quality_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        plot_files.append(quality_plot_file)
        
        print(f"   📈 Generated {len(plot_files)} visualization files")
        for plot_file in plot_files:
            print(f"      • {plot_file}")
    
    # Create comprehensive markdown report
    current_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    
    plot_section = ""
    if plot_files:
        plot_section = f"""
## 7. Visualizations

The following visualization files have been generated:

"""
        for plot_file in plot_files:
            plot_section += f"- **{plot_file}**: Analysis plots and quality metrics\n"
        
        plot_section += f"""
To view the plots, open the PNG files in any image viewer or include them in presentations.

"""

    report_content = f"""# Genomics Analysis Report

**Sample:** {sample_name}  
**Analysis Date:** {current_time}  
**Pipeline:** Flyte Genomics Pipeline v1.0

## Executive Summary

This report presents the results of quality control, read alignment, and variant calling analysis for sample **{sample_name}**.

### Overall Quality Assessment: {'✅ PASS' if alignment_metrics.mapping_rate > 80 and qc_metrics.mean_quality > 25 else '⚠️ REVIEW'}

---

## 1. Quality Control Results

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| **Total Reads** | {qc_metrics.total_reads:,} | > 10,000 | {'✅ Pass' if qc_metrics.total_reads > 10000 else '⚠️ Low'} |
| **Mean Quality Score** | {qc_metrics.mean_quality:.2f} | > 25 | {'✅ Pass' if qc_metrics.mean_quality > 25 else '⚠️ Low' if qc_metrics.mean_quality > 15 else '❌ Fail'} |
| **GC Content** | {qc_metrics.gc_content:.2f}% | 30-70% | {'✅ Pass' if 30 <= qc_metrics.gc_content <= 70 else '⚠️ Review'} |
| **Duplication Rate** | {qc_metrics.duplication_rate:.2f}% | < 30% | {'✅ Pass' if qc_metrics.duplication_rate < 30 else '⚠️ High'} |

### Quality Control Summary
- **Read Quality:** {'Excellent' if qc_metrics.mean_quality > 30 else 'Good' if qc_metrics.mean_quality > 25 else 'Poor'}
- **Data Suitability:** {'Suitable for all downstream analyses' if qc_metrics.mean_quality > 25 and 30 <= qc_metrics.gc_content <= 70 else 'May require quality filtering'}

---

## 2. Alignment Results

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| **Total Reads** | {alignment_metrics.total_reads:,} | - | - |
| **Mapped Reads** | {alignment_metrics.mapped_reads:,} | - | - |
| **Mapping Rate** | {alignment_metrics.mapping_rate:.2f}% | > 80% | {'✅ Excellent' if alignment_metrics.mapping_rate > 90 else '✅ Good' if alignment_metrics.mapping_rate > 80 else '⚠️ Low'} |
| **Mean Coverage** | {alignment_metrics.mean_coverage:.2f}x | > 20x | {'✅ Excellent' if alignment_metrics.mean_coverage > 30 else '✅ Good' if alignment_metrics.mean_coverage > 20 else '⚠️ Low'} |

### Alignment Summary
- **Mapping Efficiency:** {'Excellent alignment to reference genome' if alignment_metrics.mapping_rate > 90 else 'Good alignment efficiency' if alignment_metrics.mapping_rate > 80 else 'Suboptimal alignment - check reference genome'}
- **Coverage Assessment:** {'Sufficient depth for confident variant calling' if alignment_metrics.mean_coverage > 20 else 'Coverage may limit variant detection sensitivity'}

---

## 3. Variant Calling Results

| Variant Type | Count | Percentage | Quality Filter |
|-------------|-------|------------|----------------|
| **Total Variants** | {total_variants:,} | 100.0% | - |
| **SNPs** | {snp_count:,} | {(snp_count/total_variants*100):.1f}% | - |
| **INDELs** | {indel_count:,} | {(indel_count/total_variants*100):.1f}% | - |
| **High Quality** | {high_quality_variants:,} | {(high_quality_variants/total_variants*100):.1f}% | PASS + QUAL>30 |

### Variant Distribution
- **SNP/INDEL Ratio:** {snp_count/indel_count:.1f}:1 {'(typical)' if 3 <= snp_count/indel_count <= 10 else '(unusual - review)'}
- **Quality Rate:** {(high_quality_variants/total_variants*100):.1f}% high-quality variants
- **Variant Density:** {total_variants/alignment_metrics.mean_coverage:.1f} variants per coverage unit

---

## 4. Recommendations

### Immediate Actions
{'✅ **No action required** - All metrics within acceptable ranges' if alignment_metrics.mapping_rate > 80 and qc_metrics.mean_quality > 25 and alignment_metrics.mean_coverage > 20 else ''}
{'⚠️ **Quality Review** - Some metrics below optimal thresholds' if alignment_metrics.mapping_rate <= 80 or qc_metrics.mean_quality <= 25 or alignment_metrics.mean_coverage <= 20 else ''}

### Quality Improvements
{f'• Consider quality filtering reads with score < 20' if qc_metrics.mean_quality < 30 else '• Read quality is excellent'}
{f'• Investigate low mapping rate - check reference genome' if alignment_metrics.mapping_rate < 80 else '• Mapping rate is excellent'}
{f'• Consider increasing sequencing depth for better coverage' if alignment_metrics.mean_coverage < 25 else '• Coverage depth is sufficient'}

### Variant Analysis
• Filter variants with QUAL < 30 for high-confidence set
• Validate high-impact variants with orthogonal methods
• Consider population frequency databases for interpretation
{f'• Investigate unusual SNP/INDEL ratio' if not (3 <= snp_count/indel_count <= 10) else '• SNP/INDEL ratio is typical'}

---

## 5. Technical Details

**Pipeline Information:**
- **Software:** Flyte Genomics Pipeline
- **Version:** 1.0.0
- **Execution Mode:** {'Local simulation' if 'temp' in vcf_path else 'Production'}
- **Resource Usage:** CPU: 4 cores, Memory: 4GB total

**File Outputs:**
- Quality Control Report: Available in QC output directory
- Aligned Reads: {sample_name}_sorted.bam
- Variant Calls: {sample_name}_variants.vcf
- Analysis Report: {sample_name}_analysis_report.md

**Quality Metrics Calculation:**
- Mean Quality: Phred+33 encoding, averaged across all bases
- GC Content: Percentage of G and C nucleotides in sequences
- Mapping Rate: Percentage of reads successfully aligned to reference
- Coverage: Mean read depth across target regions

---

## 6. Conclusion

{'This sample shows excellent quality metrics suitable for all downstream genomic analyses. The data demonstrates high-quality sequencing with efficient alignment and robust variant detection.' if alignment_metrics.mapping_rate > 90 and qc_metrics.mean_quality > 30 and alignment_metrics.mean_coverage > 25 else 'This sample shows acceptable quality for genomic analysis with some areas for potential improvement noted above.'}

**Next Steps:**
1. Review variant annotations and functional impact
2. Compare variants against population databases
3. Perform downstream pathway and network analysis
4. Consider family/cohort analysis if applicable

{plot_section}---

*Report generated by Flyte Genomics Pipeline on {current_time}*
"""
    
    # Save report as markdown file
    report_file = f'{sample_name}_analysis_report.md'
    with open(report_file, 'w') as f:
        f.write(report_content)
    
    # Create summary CSV using standard library
    summary_csv = f'{sample_name}_summary.csv'
    summary_data = [
        ['metric', 'value'],
        ['total_reads', qc_metrics.total_reads],
        ['mean_quality', qc_metrics.mean_quality],
        ['gc_content', qc_metrics.gc_content],
        ['duplication_rate', qc_metrics.duplication_rate],
        ['mapping_rate', alignment_metrics.mapping_rate],
        ['mean_coverage', alignment_metrics.mean_coverage],
        ['total_variants', total_variants],
        ['snps', snp_count],
        ['indels', indel_count],
        ['high_quality_variants', high_quality_variants]
    ]
    
    with open(summary_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(summary_data)
    
    print(f"   ✅ Report Generated")
    print(f"   📄 Markdown Report: {report_file}")
    print(f"   📊 Summary CSV: {summary_csv}")
    print(f"   📈 Found {total_variants} variants ({snp_count} SNPs, {indel_count} INDELs)")
    if plot_files:
        print(f"   🎨 Generated {len(plot_files)} visualization files")
    
    return FlyteFile(path=report_file)

# Workflow Definition
@workflow
def genomics_pipeline(
    fastq_file: FlyteFile,
    reference_genome: FlyteFile, 
    sample_name: str = "sample"
) -> Tuple[FlyteDirectory, FlyteFile, FlyteFile, FlyteFile]:
    """
    Complete genomics workflow: QC -> Alignment -> Variant Calling -> Report
    
    Args:
        fastq_file: Input FASTQ sequencing file
        reference_genome: Reference genome FASTA file
        sample_name: Sample identifier for outputs
        
    Returns:
        Tuple of (QC_directory, BAM_file, VCF_file, Report_file)
    """
    
    # Step 1: Quality Control Analysis
    qc_output, qc_metrics = quality_control(fastq_file=fastq_file)
    
    # Step 2: Read Alignment to Reference
    aligned_bam, alignment_metrics = align_reads(
        fastq_file=fastq_file,
        reference_genome=reference_genome,
        sample_name=sample_name
    )
    
    # Step 3: Variant Calling
    variants_vcf = call_variants(
        bam_file=aligned_bam,
        reference_genome=reference_genome,
        sample_name=sample_name
    )
    
    # Step 4: Generate Comprehensive Report
    analysis_report = generate_report(
        qc_metrics=qc_metrics,
        alignment_metrics=alignment_metrics,
        vcf_file=variants_vcf,
        sample_name=sample_name
    )
    
    return qc_output, aligned_bam, variants_vcf, analysis_report


def create_sample_data():
    """Create sample FASTQ and reference files for testing"""
    
    print("📁 Creating sample data files...")
    
    # Create a more realistic sample FASTQ file
    fastq_content = """@SRR000001.1 HWUSI-EAS100R:6:73:941:1973 length=75
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTTTTCGTGTGGTGAAACGGT
+SRR000001.1 HWUSI-EAS100R:6:73:941:1973 length=75
CCCFFFFFHGHHGJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@SRR000001.2 HWUSI-EAS100R:6:73:941:1974 length=75
GTTTTAATGGCCGCTTAAAATTTAAATGTATTTCTGATTCTGGTTAAAAAAAGGGTTCGGCCTCCGCTCTGTAGCAGC
+SRR000001.2 HWUSI-EAS100R:6:73:941:1974 length=75
@@@DDDDDFBHFHFIIJJJJJJIIJJJJJJJJJIJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
@SRR000001.3 HWUSI-EAS100R:6:73:941:1975 length=75
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT
+SRR000001.3 HWUSI-EAS100R:6:73:941:1975 length=75
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR000001.4 HWUSI-EAS100R:6:73:941:1976 length=75
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
+SRR000001.4 HWUSI-EAS100R:6:73:941:1976 length=75
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR000001.5 HWUSI-EAS100R:6:73:941:1977 length=75
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATT
+SRR000001.5 HWUSI-EAS100R:6:73:941:1977 length=75
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR000001.6 HWUSI-EAS100R:6:73:941:1978 length=75
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
+SRR000001.6 HWUSI-EAS100R:6:73:941:1978 length=75
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SRR000001.7 HWUSI-EAS100R:6:73:941:1979 length=75
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGG
+SRR000001.7 HWUSI-EAS100R:6:73:941:1979 length=75
CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
"""
    
    # Create a sample reference genome file
    reference_content = """>chr1 chromosome 1, complete sequence
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
>chr2 chromosome 2, complete sequence
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT
TTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
ACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACAC
TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA
"""
    
    # Write files to current directory
    with open("sample.fastq", "w") as f:
        f.write(fastq_content)
    
    with open("reference.fa", "w") as f:
        f.write(reference_content)
    
    print("   ✅ sample.fastq created (7 reads, 75bp each)")
    print("   ✅ reference.fa created (2 chromosomes)")
    
    return "sample.fastq", "reference.fa"


def main():
    """Main function with usage examples"""
    
    # Additional warning suppression for runtime
    import warnings
    warnings.filterwarnings("ignore", category=UserWarning, module="click")
    
    print("🧬 Flyte Bioinformatics Workflow")
    print("=" * 50)
    
    # Create sample data if files don't exist
    if not os.path.exists("sample.fastq") or not os.path.exists("reference.fa"):
        fastq_file, ref_file = create_sample_data()
    else:
        fastq_file, ref_file = "sample.fastq", "reference.fa"
        print("📁 Using existing sample files")
    
    print(f"\n📋 Usage Examples:")
    print("=" * 50)
    
    print("\n1️⃣ Run workflow locally (recommended for testing):")
    print(f"   pyflyte run {__file__} genomics_pipeline \\")
    print(f"       --fastq_file {fastq_file} \\")
    print(f"       --reference_genome {ref_file} \\")
    print(f"       --sample_name 'test_sample'")
    
    print("\n2️⃣ Run individual tasks:")
    print(f"   pyflyte run {__file__} quality_control --fastq_file {fastq_file}")
    print(f"   pyflyte run {__file__} align_reads --fastq_file {fastq_file} --reference_genome {ref_file}")
    
    print("\n3️⃣ Suppress CLI warnings (if needed):")
    print(f"   PYTHONWARNINGS=ignore::UserWarning pyflyte run {__file__} genomics_pipeline \\")
    print(f"       --fastq_file {fastq_file} --reference_genome {ref_file}")
    
    print("\n4️⃣ Run with Flyte cluster (auto-installs dependencies):")
    print("   # First start the cluster:")
    print("   flytectl demo start")
    print("")
    print("   # Register workflow:")
    print(f"   pyflyte --config ~/.flyte/config-sandbox.yaml register {__file__}")
    print("")
    print("   # Run on cluster:")
    print(f"   pyflyte --config ~/.flyte/config-sandbox.yaml run --remote {__file__} genomics_pipeline \\")
    print(f"       --fastq_file {fastq_file} --reference_genome {ref_file}")
    
    print("\n5️⃣ View Flyte UI (if using sandbox):")
    print("   Open: http://localhost:30080")
    
    print("\n📦 Requirements:")
    print("   pip install \"flytekit[all]\"")
    print("   Optional: matplotlib (for visualizations)")
    print("   Optional: flytectl (for sandbox)")
    print("   Optional: Docker (for cluster execution)")
    
    print("\n📄 Expected Outputs:")
    print("   • QC report directory with HTML analysis")
    print("   • Aligned BAM file (simulated)")
    print("   • VCF file with variant calls")
    print("   • Comprehensive markdown analysis report")
    print("   • Summary CSV with key metrics")
    print("   • Visualization plots (PNG files, if matplotlib available)")
    
    print("\n💡 Cluster Execution Notes:")
    print("   • Local execution uses your Python environment")
    print("   • Cluster execution uses only standard Python libraries")
    print("   • Matplotlib auto-installed if needed for visualizations")
    print("   • Compatible with any Flyte sandbox or cluster")
    
    print("\n🐛 Troubleshooting:")
    print("   • No more ModuleNotFoundError issues!")
    print("   • File not found in cluster? Use absolute paths for inputs")
    print("   • CLI warnings? Use PYTHONWARNINGS=ignore environment variable")


if __name__ == "__main__":
    # Set environment variable to suppress warnings
    os.environ.setdefault('PYTHONWARNINGS', 'ignore::UserWarning')
    main()

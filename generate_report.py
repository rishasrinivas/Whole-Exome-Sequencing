

import argparse
import sys
import os
import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from jinja2 import Template
import pysam
import re
from collections import defaultdict, Counter


class VCFAnnotationParser:
    """Parser for VEP-annotated VCF files"""
    
    def __init__(self, vcf_path):
        self.vcf_path = vcf_path
        self.variants = []
        self.sample_name = None
        self.vep_header = {}
        
    def parse_vcf(self):
        """Parse the VCF file and extract variant information"""
        try:
            vcf = pysam.VariantFile(self.vcf_path)
            
            # Get sample name
            self.sample_name = list(vcf.header.samples)[0] if vcf.header.samples else "Unknown"
            
            # Parse VEP header info
            self._parse_vep_header(vcf.header)
            
            # Parse variants
            for record in vcf:
                variant_info = self._parse_variant_record(record)
                if variant_info:
                    self.variants.append(variant_info)
                    
        except Exception as e:
            print(f"Error parsing VCF file: {e}")
            sys.exit(1)
            
        return self.variants
    
    def _parse_vep_header(self, header):
        """Parse VEP header information"""
        for rec in header.records:
            if rec.key == 'INFO' and rec.id == 'CSQ':
                # Extract CSQ field format
                desc = rec.description
                format_match = re.search(r'Format: (.+)', desc)
                if format_match:
                    fields = [f.strip() for f in format_match.group(1).split('|')]
                    self.vep_header = {i: field for i, field in enumerate(fields)}
    
    def _parse_variant_record(self, record):
        """Parse individual variant record"""
        try:
            # Basic variant information
            variant = {
                'chrom': record.chrom,
                'pos': record.pos,
                'id': record.id if record.id else '.',
                'ref': record.ref,
                'alt': ','.join([str(alt) for alt in record.alts]),
                'qual': record.qual if record.qual else 0,
                'filter': ','.join(record.filter) if record.filter else 'PASS',
            }
            
            # Sample genotype information
            if record.samples:
                sample = list(record.samples.values())[0]
                variant.update({
                    'gt': sample.get('GT', './.'),
                    'dp': sample.get('DP', 0),
                    'gq': sample.get('GQ', 0),
                    'ad': sample.get('AD', [0, 0]),
                    'vaf': self._calculate_vaf(sample.get('AD', [0, 0]))
                })
            
            # Parse VEP annotations
            vep_annotations = self._parse_vep_annotation(record)
            variant['annotations'] = vep_annotations
            
            # Determine most severe consequence
            variant['most_severe_consequence'] = self._get_most_severe_consequence(vep_annotations)
            
            # Extract clinical significance
            variant['clinical_significance'] = self._extract_clinical_significance(vep_annotations)
            
            return variant
            
        except Exception as e:
            print(f"Error parsing variant record: {e}")
            return None
    
    def _calculate_vaf(self, ad):
        """Calculate variant allele frequency"""
        if len(ad) >= 2 and sum(ad) > 0:
            return ad[1] / sum(ad)
        return 0
    
    def _parse_vep_annotation(self, record):
        """Parse VEP CSQ annotations"""
        annotations = []
        
        if 'CSQ' in record.info:
            for csq in record.info['CSQ']:
                csq_fields = csq.split('|')
                annotation = {}
                
                for i, value in enumerate(csq_fields):
                    if i in self.vep_header:
                        field_name = self.vep_header[i]
                        annotation[field_name] = value if value else None
                
                annotations.append(annotation)
        
        return annotations
    
    def _get_most_severe_consequence(self, annotations):
        """Determine the most severe consequence from VEP annotations"""
        severity_order = [
            'transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
            'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost',
            'transcript_amplification', 'inframe_insertion', 'inframe_deletion',
            'missense_variant', 'protein_altering_variant', 'splice_region_variant',
            'incomplete_terminal_codon_variant', 'start_retained_variant',
            'stop_retained_variant', 'synonymous_variant', 'coding_sequence_variant',
            'mature_miRNA_variant', '5_prime_UTR_variant', '3_prime_UTR_variant',
            'non_coding_transcript_exon_variant', 'intron_variant',
            'NMD_transcript_variant', 'non_coding_transcript_variant',
            'upstream_gene_variant', 'downstream_gene_variant', 'TFBS_ablation',
            'TFBS_amplification', 'TF_binding_site_variant', 'regulatory_region_ablation',
            'regulatory_region_amplification', 'feature_elongation',
            'regulatory_region_variant', 'feature_truncation', 'intergenic_variant'
        ]
        
        most_severe = 'intergenic_variant'
        
        for annotation in annotations:
            consequence = annotation.get('Consequence', '')
            if consequence in severity_order:
                if severity_order.index(consequence) < severity_order.index(most_severe):
                    most_severe = consequence
        
        return most_severe
    
    def _extract_clinical_significance(self, annotations):
        """Extract clinical significance from ClinVar annotations"""
        clinical_significances = []
        
        for annotation in annotations:
            clin_sig = annotation.get('CLIN_SIG', '')
            if clin_sig:
                clinical_significances.append(clin_sig)
        
        # Return the most significant classification
        priority = ['pathogenic', 'likely_pathogenic', 'uncertain_significance', 
                   'likely_benign', 'benign']
        
        for sig in priority:
            if any(sig.replace('_', ' ') in cs.lower() for cs in clinical_significances):
                return sig.replace('_', ' ').title()
        
        return 'Unknown' if clinical_significances else None


class ClinicalReportGenerator:
    """Generate clinical reports from parsed VCF data"""
    
    def __init__(self, variants, sample_name, output_dir):
        self.variants = variants
        self.sample_name = sample_name
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Filter variants for clinical relevance
        self.pathogenic_variants = self._filter_pathogenic_variants()
        self.variant_summary = self._generate_variant_summary()
        
    def _filter_pathogenic_variants(self):
        """Filter variants with high clinical significance"""
        pathogenic = []
        
        for variant in self.variants:
            clin_sig = variant.get('clinical_significance', '')
            if clin_sig and ('Pathogenic' in clin_sig or 'Likely Pathogenic' in clin_sig):
                pathogenic.append(variant)
        
        return sorted(pathogenic, key=lambda x: (x['chrom'], x['pos']))
    
    def _generate_variant_summary(self):
        """Generate summary statistics"""
        summary = {
            'total_variants': len(self.variants),
            'pathogenic_variants': len([v for v in self.variants 
                                      if v.get('clinical_significance') == 'Pathogenic']),
            'likely_pathogenic_variants': len([v for v in self.variants 
                                             if v.get('clinical_significance') == 'Likely Pathogenic']),
            'uncertain_significance': len([v for v in self.variants 
                                         if v.get('clinical_significance') == 'Uncertain Significance']),
            'benign_variants': len([v for v in self.variants 
                                  if v.get('clinical_significance') in ['Benign', 'Likely Benign']]),
        }
        
        # Consequence summary
        consequences = [v['most_severe_consequence'] for v in self.variants]
        summary['consequence_counts'] = Counter(consequences)
        
        return summary
    
    def generate_html_report(self):
        """Generate HTML clinical report"""
        html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>WES Clinical Report - {{ sample_name }}</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; line-height: 1.6; }
        .header { background: #f4f4f4; padding: 20px; border-radius: 5px; margin-bottom: 20px; }
        .section { margin-bottom: 30px; }
        .pathogenic { background: #ffebee; border-left: 4px solid #f44336; padding: 10px; }
        .likely-pathogenic { background: #fff3e0; border-left: 4px solid #ff9800; padding: 10px; }
        table { border-collapse: collapse; width: 100%; margin-top: 10px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #f2f2f2; }
        .high-impact { color: #d32f2f; font-weight: bold; }
        .moderate-impact { color: #f57c00; }
        .summary-box { background: #e3f2fd; padding: 15px; border-radius: 5px; margin: 10px 0; }
        .gene { font-weight: bold; color: #1976d2; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Whole Exome Sequencing Clinical Report</h1>
        <p><strong>Sample:</strong> {{ sample_name }}</p>
        <p><strong>Report Generated:</strong> {{ report_date }}</p>
        <p><strong>Pipeline Version:</strong> 1.0.0</p>
    </div>

    <div class="section">
        <h2>Executive Summary</h2>
        <div class="summary-box">
            <p><strong>Total Variants:</strong> {{ summary.total_variants }}</p>
            <p><strong>Pathogenic Variants:</strong> {{ summary.pathogenic_variants }}</p>
            <p><strong>Likely Pathogenic Variants:</strong> {{ summary.likely_pathogenic_variants }}</p>
            <p><strong>Variants of Uncertain Significance:</strong> {{ summary.uncertain_significance }}</p>
            <p><strong>Benign/Likely Benign Variants:</strong> {{ summary.benign_variants }}</p>
        </div>
    </div>

    {% if pathogenic_variants %}
    <div class="section">
        <h2>Clinically Significant Variants</h2>
        <p>The following variants have been classified as Pathogenic or Likely Pathogenic and require clinical attention:</p>
        
        {% for variant in pathogenic_variants %}
        <div class="{% if variant.clinical_significance == 'Pathogenic' %}pathogenic{% else %}likely-pathogenic{% endif %}">
            <h3>Variant {{ loop.index }}</h3>
            <table>
                <tr><th>Location</th><td>{{ variant.chrom }}:{{ variant.pos }}</td></tr>
                <tr><th>Change</th><td>{{ variant.ref }} &rarr; {{ variant.alt }}</td></tr>
                <tr><th>Clinical Significance</th><td class="high-impact">{{ variant.clinical_significance }}</td></tr>
                <tr><th>Consequence</th><td>{{ variant.most_severe_consequence.replace('_', ' ').title() }}</td></tr>
                <tr><th>Genotype</th><td>{{ variant.gt }}</td></tr>
                <tr><th>Depth</th><td>{{ variant.dp }}</td></tr>
                <tr><th>Quality</th><td>{{ variant.gq }}</td></tr>
                <tr><th>VAF</th><td>{{ "%.2f"|format(variant.vaf * 100) }}%</td></tr>
            </table>
            
            {% if variant.annotations %}
            <h4>Gene Information</h4>
            <table>
                <tr><th>Gene</th><th>Transcript</th><th>HGVSc</th><th>HGVSp</th><th>Impact</th></tr>
                {% for annotation in variant.annotations[:3] %}
                <tr>
                    <td class="gene">{{ annotation.SYMBOL or 'N/A' }}</td>
                    <td>{{ annotation.Feature or 'N/A' }}</td>
                    <td>{{ annotation.HGVSc or 'N/A' }}</td>
                    <td>{{ annotation.HGVSp or 'N/A' }}</td>
                    <td>{{ annotation.IMPACT or 'N/A' }}</td>
                </tr>
                {% endfor %}
            </table>
            {% endif %}
        </div>
        {% endfor %}
    </div>
    {% else %}
    <div class="section">
        <h2>Clinically Significant Variants</h2>
        <p>No Pathogenic or Likely Pathogenic variants were found in this sample.</p>
    </div>
    {% endif %}

    <div class="section">
        <h2>Appendix: Full Variant Summary</h2>
        <p>This section provides a complete breakdown of all identified variants.</p>
        <div class="summary-box">
            <h4>Variant Consequence Counts</h4>
            <pre>{{ summary.consequence_counts | tojson(indent=2) }}</pre>
        </div>
    </div>
    
    <footer>
        <p>Generated by WES Pipeline Team</p>
    </footer>

</body>
</html>
"""
        template = Template(html_template)
        report_date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        html_output = template.render(
            sample_name=self.sample_name,
            report_date=report_date,
            summary=self.variant_summary,
            pathogenic_variants=self.pathogenic_variants,
        )
        
        output_file = self.output_dir / f"{self.sample_name}_clinical_report.html"
        with open(output_file, 'w') as f:
            f.write(html_output)
        
        print(f"Clinical report generated at: {output_file}")


def main():
    """Main function to parse arguments and run the report generator"""
    parser = argparse.ArgumentParser(
        description="Generate a clinical report from a VEP-annotated VCF file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("vcf_file", help="Path to the VEP-annotated VCF file.")
    parser.add_argument(
        "-o", "--output_dir", 
        default="clinical_report", 
        help="Output directory for the generated report (default: 'clinical_report')."
    )

    args = parser.parse_args()

    # Check if VCF file exists
    if not Path(args.vcf_file).is_file():
        print(f"Error: VCF file not found at '{args.vcf_file}'")
        sys.exit(1)

    print("Parsing VCF file...")
    parser = VCFAnnotationParser(args.vcf_file)
    variants = parser.parse_vcf()

    if not variants:
        print("No variants found in the VCF file. Exiting.")
        sys.exit(0)
    
    print(f"Found {len(variants)} variants for sample '{parser.sample_name}'.")

    print("Generating clinical report...")
    report_generator = ClinicalReportGenerator(
        variants=variants,
        sample_name=parser.sample_name,
        output_dir=args.output_dir
    )
    report_generator.generate_html_report()
    print("Report generation complete.")


if __name__ == "__main__":
    main()

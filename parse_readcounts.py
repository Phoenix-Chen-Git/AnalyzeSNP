
# def parse_snp_info(snp_entries):
#     counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
#     for entry in snp_entries:
#         if entry:  # Ensure entry is not empty
#             parts = entry.split(':')
#             base = parts[0]
#             if base in counts:  # Only process A, T, G, C
#                 reads = int(parts[1])
#                 counts[base] = reads
#     total_reads = sum(counts.values())
#     return counts, total_reads

# def calculate_mut_freq(ref_base, counts, total_reads):
#     if total_reads == 0:
#         return 0
#     ref_reads = counts.get(ref_base, 0)
#     mut_reads = total_reads - ref_reads
#     return mut_reads / total_reads if total_reads > 0 else 0

# def parse_readcounts(input_file, output_prefix):
#     txt_file = f"{output_prefix}.txt"
#     csv_file = f"{output_prefix}.csv"

#     with open(input_file, 'r') as f_in, \
#          open(txt_file, 'w') as f_txt, \
#          open(csv_file, 'w') as f_csv:
        
#         # 更新CSV头，添加12列突变列
#         header = "chrom,position,ref_base,A,T,G,C,ref>A,ref>T,ref>G,ref>C,mut_freq," + \
#                  "A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G"
#         f_txt.write('\t'.join(header.split(',')) + '\n')
#         f_csv.write(header + '\n')

#         f_in.readline()  # 跳过输入头

#         for line in f_in:
#             fields = line.strip().split('\t')
#             chrom = fields[0]
#             position = fields[1]
#             ref_base = fields[2].upper()
#             snp_entries = fields[5:]  # All columns starting from the 6th
#             counts, total_reads = parse_snp_info(snp_entries)
            
#             # 计算现有频率
#             if total_reads > 0:
#                 ref_to_a = counts['A'] / total_reads if ref_base != 'A' else 0
#                 ref_to_t = counts['T'] / total_reads if ref_base != 'T' else 0
#                 ref_to_g = counts['G'] / total_reads if ref_base != 'G' else 0
#                 ref_to_c = counts['C'] / total_reads if ref_base != 'C' else 0
#             else:
#                 ref_to_a = ref_to_t = ref_to_g = ref_to_c = 0
#             mut_freq = calculate_mut_freq(ref_base, counts, total_reads)

#             # 新增12列突变频率（基于ref_base决定哪些列有值）
#             mutation_cols = {
#                 'A>T': counts['T'] / total_reads if ref_base == 'A' and total_reads > 0 else 0,
#                 'A>G': counts['G'] / total_reads if ref_base == 'A' and total_reads > 0 else 0,
#                 'A>C': counts['C'] / total_reads if ref_base == 'A' and total_reads > 0 else 0,
#                 'T>A': counts['A'] / total_reads if ref_base == 'T' and total_reads > 0 else 0,
#                 'T>G': counts['G'] / total_reads if ref_base == 'T' and total_reads > 0 else 0,
#                 'T>C': counts['C'] / total_reads if ref_base == 'T' and total_reads > 0 else 0,
#                 'G>A': counts['A'] / total_reads if ref_base == 'G' and total_reads > 0 else 0,
#                 'G>T': counts['T'] / total_reads if ref_base == 'G' and total_reads > 0 else 0,
#                 'G>C': counts['C'] / total_reads if ref_base == 'G' and total_reads > 0 else 0,
#                 'C>A': counts['A'] / total_reads if ref_base == 'C' and total_reads > 0 else 0,
#                 'C>T': counts['T'] / total_reads if ref_base == 'C' and total_reads > 0 else 0,
#                 'C>G': counts['G'] / total_reads if ref_base == 'C' and total_reads > 0 else 0
#             }

#             # 构建CSV行
#             csv_row = [
#                 chrom, position, ref_base,
#                 str(counts['A']), str(counts['T']),
#                 str(counts['G']), str(counts['C']),
#                 f"{ref_to_a:.6f}", f"{ref_to_t:.6f}",
#                 f"{ref_to_g:.6f}", f"{ref_to_c:.6f}",
#                 f"{mut_freq:.6f}"
#             ] + [f"{mutation_cols[col]:.6f}" for col in mutation_cols]
#             f_csv.write(','.join(csv_row) + '\n')

#             # 构建TXT行
#             txt_row = [
#                 chrom, position, ref_base,
#                 str(counts['A']), str(counts['T']),
#                 str(counts['G']), str(counts['C']),
#                 f"{ref_to_a:.6f}", f"{ref_to_t:.6f}",
#                 f"{ref_to_g:.6f}", f"{ref_to_c:.6f}",
#                 f"{mut_freq:.6f}"
#             ] + [f"{mutation_cols[col]:.6f}" for col in mutation_cols]
#             f_txt.write('\t'.join(txt_row) + '\n')

# def generate_mutation_matrix(sample_name, output_dir):
#     """生成4x4碱基替换频率矩阵（纯小数格式）"""
#     input_csv = f"{output_dir}/{sample_name}_full_base.csv"
#     output_file = f"{output_dir}/{sample_name}_mutation_matrix.csv"

#     # 初始化数据结构 (ref_base -> sub_base)
#     bases = ['A', 'T', 'G', 'C']
#     transition_data = {
#         ref: {sub: {'total': 0.0, 'weight': 0} for sub in bases if sub != ref}
#         for ref in bases
#     }

#     # 更新列号映射表（基于新CSV结构）
#     column_map = {
#         'A': [('T', 12), ('G', 13), ('C', 14)],  # A>T(12),A>G(13),A>C(14)
#         'T': [('A', 15), ('G', 16), ('C', 17)],  # T>A(15),T>G(16),T>C(17)
#         'G': [('A', 18), ('T', 19), ('C', 20)],  # G>A(18),G>T(19),G>C(20)
#         'C': [('A', 21), ('T', 22), ('G', 23)]   # C>A(21),C>T(22),C>G(23)
#     }

#     # 读取CSV文件
#     with open(input_csv, 'r') as f:
#         next(f)  # 跳过标题行
#         for line in f:
#             row = line.strip().split(',')
#             ref_base = row[2].upper()
            
#             # 跳过无效碱基
#             if ref_base not in column_map:
#                 continue

#             # 计算总reads数
#             try:
#                 counts = [int(row[i]) for i in [3,4,5,6]]  # A(3),T(4),G(5),C(6)
#                 total = sum(counts)
#                 if total == 0:
#                     continue
#             except ValueError:
#                 continue

#             # 累加加权值
#             for sub_base, col_idx in column_map[ref_base]:
#                 freq = float(row[col_idx])
#                 transition_data[ref_base][sub_base]['total'] += freq * total
#                 transition_data[ref_base][sub_base]['weight'] += total

#     # 计算加权平均
#     mutation_matrix = {}
#     for ref in bases:
#         mutation_matrix[ref] = {}
#         for sub in bases:
#             if ref == sub:
#                 mutation_matrix[ref][sub] = 0.0
#             else:
#                 data = transition_data[ref].get(sub, {'total':0, 'weight':0})
#                 if data['weight'] > 0:
#                     mutation_matrix[ref][sub] = data['total'] / data['weight']
#                 else:
#                     mutation_matrix[ref][sub] = 0.0

#     # 写入CSV文件
#     with open(output_file, 'w') as f:
#         # 写入表头
#         f.write("Reference\\Substitute,A,T,G,C\n")
        
#         # 写入数据行
#         for ref in bases:
#             row = [ref]
#             for sub in bases:
#                 value = mutation_matrix[ref][sub]
#                 formatted = "0" if value == 0 else f"{value:.6f}"
#                 row.append(formatted)
#             f.write(",".join(row) + "\n")

# if __name__ == "__main__":
#     input_file = sys.argv[1]
#     sample_name = sys.argv[2]
#     output_dir = sys.argv[3]
    
#     output_prefix = f"{output_dir}/{sample_name}_full_base"
#     parse_readcounts(input_file, output_prefix)
#     generate_mutation_matrix(sample_name, output_dir)

import sys
import math
from statistics import median

def parse_snp_info(snp_entries):
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for entry in snp_entries:
        if entry:  # Ensure entry is not empty
            parts = entry.split(':')
            base = parts[0]
            if base in counts:  # Only process A, T, G, C
                reads = int(parts[1])
                counts[base] = reads
    total_reads = sum(counts.values())
    return counts, total_reads

def calculate_mut_freq(ref_base, counts, total_reads):
    if total_reads == 0:
        return 0
    ref_reads = counts.get(ref_base, 0)
    mut_reads = total_reads - ref_reads
    return mut_reads / total_reads if total_reads > 0 else 0

def parse_readcounts(input_file, output_prefix):
    txt_file = f"{output_prefix}.txt"
    csv_file = f"{output_prefix}.csv"

    with open(input_file, 'r') as f_in, \
         open(txt_file, 'w') as f_txt, \
         open(csv_file, 'w') as f_csv:
        
        # 更新CSV头，添加12列突变列
        header = "chrom,position,ref_base,A,T,G,C,ref>A,ref>T,ref>G,ref>C,mut_freq," + \
                 "A>T,A>G,A>C,T>A,T>G,T>C,G>A,G>T,G>C,C>A,C>T,C>G"
        f_txt.write('\t'.join(header.split(',')) + '\n')
        f_csv.write(header + '\n')

        f_in.readline()  # 跳过输入头

        for line in f_in:
            fields = line.strip().split('\t')
            chrom = fields[0]
            position = fields[1]
            ref_base = fields[2].upper()
            snp_entries = fields[5:]  # All columns starting from the 6th
            counts, total_reads = parse_snp_info(snp_entries)
            
            # 计算现有频率
            if total_reads > 0:
                ref_to_a = counts['A'] / total_reads if ref_base != 'A' else 0
                ref_to_t = counts['T'] / total_reads if ref_base != 'T' else 0
                ref_to_g = counts['G'] / total_reads if ref_base != 'G' else 0
                ref_to_c = counts['C'] / total_reads if ref_base != 'C' else 0
            else:
                ref_to_a = ref_to_t = ref_to_g = ref_to_c = 0
            mut_freq = calculate_mut_freq(ref_base, counts, total_reads)

            # 新增12列突变频率（基于ref_base决定哪些列有值）
            mutation_cols = {
                'A>T': counts['T'] / total_reads if ref_base == 'A' and total_reads > 0 else 0,
                'A>G': counts['G'] / total_reads if ref_base == 'A' and total_reads > 0 else 0,
                'A>C': counts['C'] / total_reads if ref_base == 'A' and total_reads > 0 else 0,
                'T>A': counts['A'] / total_reads if ref_base == 'T' and total_reads > 0 else 0,
                'T>G': counts['G'] / total_reads if ref_base == 'T' and total_reads > 0 else 0,
                'T>C': counts['C'] / total_reads if ref_base == 'T' and total_reads > 0 else 0,
                'G>A': counts['A'] / total_reads if ref_base == 'G' and total_reads > 0 else 0,
                'G>T': counts['T'] / total_reads if ref_base == 'G' and total_reads > 0 else 0,
                'G>C': counts['C'] / total_reads if ref_base == 'G' and total_reads > 0 else 0,
                'C>A': counts['A'] / total_reads if ref_base == 'C' and total_reads > 0 else 0,
                'C>T': counts['T'] / total_reads if ref_base == 'C' and total_reads > 0 else 0,
                'C>G': counts['G'] / total_reads if ref_base == 'C' and total_reads > 0 else 0
            }

            # 构建CSV行
            csv_row = [
                chrom, position, ref_base,
                str(counts['A']), str(counts['T']),
                str(counts['G']), str(counts['C']),
                f"{ref_to_a:.6f}", f"{ref_to_t:.6f}",
                f"{ref_to_g:.6f}", f"{ref_to_c:.6f}",
                f"{mut_freq:.6f}"
            ] + [f"{mutation_cols[col]:.6f}" for col in mutation_cols]
            f_csv.write(','.join(csv_row) + '\n')

            # 构建TXT行
            txt_row = [
                chrom, position, ref_base,
                str(counts['A']), str(counts['T']),
                str(counts['G']), str(counts['C']),
                f"{ref_to_a:.6f}", f"{ref_to_t:.6f}",
                f"{ref_to_g:.6f}", f"{ref_to_c:.6f}",
                f"{mut_freq:.6f}"
            ] + [f"{mutation_cols[col]:.6f}" for col in mutation_cols]
            f_txt.write('\t'.join(txt_row) + '\n')

def generate_mutation_matrix(sample_name, output_dir):
    """生成4x4碱基替换频率矩阵（算术平均、几何平均、中位数）"""
    input_csv = f"{output_dir}/{sample_name}_full_base.csv"
    output_file_mean = f"{output_dir}/{sample_name}_mutation_matrix.csv"
    output_file_geometric = f"{output_dir}/{sample_name}_mutation_matrix_geometric.csv"
    output_file_median = f"{output_dir}/{sample_name}_mutation_matrix_median.csv"

    # 初始化数据结构
    bases = ['A', 'T', 'G', 'C']
    transition_data = {
        ref: {sub: {'total': 0.0, 'weight': 0, 'freqs': []} for sub in bases if sub != ref}
        for ref in bases
    }

    # 更新列号映射表（基于CSV结构）
    column_map = {
        'A': [('T', 12), ('G', 13), ('C', 14)],  # A>T(12),A>G(13),A>C(14)
        'T': [('A', 15), ('G', 16), ('C', 17)],  # T>A(15),T>G(16),T>C(17)
        'G': [('A', 18), ('T', 19), ('C', 20)],  # G>A(18),G>T(19),G>C(20)
        'C': [('A', 21), ('T', 22), ('G', 23)]   # C>A(21),C>T(22),C>G(23)
    }

    # 读取CSV文件
    with open(input_csv, 'r') as f:
        next(f)  # 跳过标题行
        for line in f:
            row = line.strip().split(',')
            ref_base = row[2].upper()
            
            # 跳过无效碱基
            if ref_base not in column_map:
                continue

            # 计算总reads数
            try:
                counts = [int(row[i]) for i in [3,4,5,6]]  # A(3),T(4),G(5),C(6)
                total = sum(counts)
                if total == 0:
                    continue
            except ValueError:
                continue

            # 累加加权值并收集非零频率
            for sub_base, col_idx in column_map[ref_base]:
                freq = float(row[col_idx])
                transition_data[ref_base][sub_base]['total'] += freq * total
                transition_data[ref_base][sub_base]['weight'] += total
                if freq > 0:
                    transition_data[ref_base][sub_base]['freqs'].append(freq)

    # 计算算术平均
    mutation_matrix_mean = {}
    for ref in bases:
        mutation_matrix_mean[ref] = {}
        for sub in bases:
            if ref == sub:
                mutation_matrix_mean[ref][sub] = 0.0
            else:
                data = transition_data[ref].get(sub, {'total':0, 'weight':0})
                if data['weight'] > 0:
                    mutation_matrix_mean[ref][sub] = data['total'] / data['weight']
                else:
                    mutation_matrix_mean[ref][sub] = 0.0

    # 计算几何平均
    mutation_matrix_geometric = {}
    for ref in bases:
        mutation_matrix_geometric[ref] = {}
        for sub in bases:
            if ref == sub:
                mutation_matrix_geometric[ref][sub] = 0.0
            else:
                freqs = transition_data[ref].get(sub, {'freqs': []})['freqs']
                if freqs:
                    # 几何平均 = exp(平均(ln(频率)))
                    log_mean = sum(math.log(f) for f in freqs) / len(freqs)
                    mutation_matrix_geometric[ref][sub] = math.exp(log_mean)
                else:
                    mutation_matrix_geometric[ref][sub] = 0.0

    # 计算中位数
    mutation_matrix_median = {}
    for ref in bases:
        mutation_matrix_median[ref] = {}
        for sub in bases:
            if ref == sub:
                mutation_matrix_median[ref][sub] = 0.0
            else:
                freqs = transition_data[ref].get(sub, {'freqs': []})['freqs']
                if freqs:
                    mutation_matrix_median[ref][sub] = median(freqs)
                else:
                    mutation_matrix_median[ref][sub] = 0.0

    # 写入算术平均矩阵
    with open(output_file_mean, 'w') as f:
        f.write("Reference\\Substitute,A,T,G,C\n")
        for ref in bases:
            row = [ref]
            for sub in bases:
                value = mutation_matrix_mean[ref][sub]
                formatted = "0" if value == 0 else f"{value:.6f}"
                row.append(formatted)
            f.write(",".join(row) + "\n")

    # 写入几何平均矩阵
    with open(output_file_geometric, 'w') as f:
        f.write("Reference\\Substitute,A,T,G,C\n")
        for ref in bases:
            row = [ref]
            for sub in bases:
                value = mutation_matrix_geometric[ref][sub]
                formatted = "0" if value == 0 else f"{value:.6f}"
                row.append(formatted)
            f.write(",".join(row) + "\n")

    # 写入中位数矩阵
    with open(output_file_median, 'w') as f:
        f.write("Reference\\Substitute,A,T,G,C\n")
        for ref in bases:
            row = [ref]
            for sub in bases:
                value = mutation_matrix_median[ref][sub]
                formatted = "0" if value == 0 else f"{value:.6f}"
                row.append(formatted)
            f.write(",".join(row) + "\n")

def extract_nonzero_mutation(sample_name, output_dir):
    """Extract non-zero mutation frequencies for all mutation types, aligned by columns."""
    input_csv = f"{output_dir}/{sample_name}_full_base.csv"
    output_file = f"{output_dir}/{sample_name}_all_mutations_nonzero.txt"
    
    # Define mutation types and their column indices
    mutation_columns = {
        'A>T': 12, 'A>G': 13, 'A>C': 14,
        'T>A': 15, 'T>G': 16, 'T>C': 17,
        'G>A': 18, 'G>T': 19, 'G>C': 20,
        'C>A': 21, 'C>T': 22, 'C>G': 23
    }
    
    # Collect non-zero frequencies for each mutation type
    frequencies = {mut: [] for mut in mutation_columns}
    with open(input_csv, 'r') as f:
        next(f)  # Skip header
        for line in f:
            row = line.strip().split(',')
            for mutation, col_idx in mutation_columns.items():
                try:
                    freq = float(row[col_idx])
                    if freq > 0:
                        frequencies[mutation].append(f"{freq:.6f}")
                except (ValueError, IndexError):
                    continue
    
    # Determine the maximum number of frequencies for any mutation type
    max_rows = max(len(freqs) for freqs in frequencies.values())
    
    # Write to output file with header and aligned columns
    with open(output_file, 'w') as f:
        # Write header
        header = ','.join(mutation_columns.keys())
        f.write(header + '\n')
        
        # Write data rows
        for i in range(max_rows):
            row = []
            for mutation in mutation_columns:
                # Get the frequency at index i, or empty string if none
                freq_list = frequencies[mutation]
                freq = freq_list[i] if i < len(freq_list) else ""
                row.append(freq)
            f.write(','.join(row) + '\n')

if __name__ == "__main__":
    input_file = sys.argv[1]
    sample_name = sys.argv[2]
    output_dir = sys.argv[3]
    
    output_prefix = f"{output_dir}/{sample_name}_full_base"
    parse_readcounts(input_file, output_prefix)
    generate_mutation_matrix(sample_name, output_dir)
    extract_nonzero_mutation(sample_name, output_dir)
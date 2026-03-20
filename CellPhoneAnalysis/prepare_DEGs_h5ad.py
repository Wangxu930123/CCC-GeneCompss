"""
DEGs文件生成与CellPhoneDB分析脚本

功能:
1. 从h5ad数据自动生成符合CellPhoneDB标准的DEGs文件
2. 运行CellPhoneDB细胞互作分析
3. 支持多种差异分析方法（Wilcoxon/t-test/logreg）
4. 智能检测基因标识符类型（Ensembl/HGNC/基因名）

DEGs文件格式:
- CellPhoneDB标准格式: gene, cell_type1, cell_type2, ...
  每行一个基因，每列是该细胞类型的差异基因列表（逗号分隔，按log2FC降序）

使用方法:
    # 仅生成DEGs文件
    python prepare_DEGs_h5ad.py --mode degs_only --h5ad data.h5ad --outdir ./output

    # 完整流程，自动生成DEGs
    python prepare_DEGs_h5ad.py --mode full --h5ad data.h5ad --generate_degs --outdir ./output

    # 完整流程，使用已有DEGs
    python prepare_DEGs_h5ad.py --mode full --h5ad data.h5ad --degs DEGs.tsv --microenv microenv.tsv

作者: GeneCompass团队
"""

import os
import logging
import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import scipy
from cellphonedb.src.core.methods import cpdb_degs_analysis_method
import multiprocessing

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def generate_degs_file(adata, output_path, groupby='cell_type', method='wilcoxon', n_top_genes=250):
    """
    从h5ad数据生成符合CellPhoneDB标准的DEGs文件

    参数:
        adata: AnnData对象
        output_path: 输出DEGs文件路径
        groupby: 分组列名
        method: 差异分析方法 ('wilcoxon', 't-test', 'logreg')
        n_top_genes: 每个细胞类型选取的差异基因数量

    CellPhoneDB DEGs文件格式:
        - 第1列: 基因名称 (所有行相同)
        - 第2列+: 每个细胞类型的差异基因列表 (按log2FC降序排列)
    """
    try:
        logger.info(f"开始生成DEGs文件，方法: {method}")

        # 1. 标准化数据
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # 2. 计算差异表达基因
        logger.info(f"计算差异表达基因 ({method})...")
        sc.tl.rank_genes_groups(
            adata,
            groupby=groupby,
            method=method,
            corr_method='bonferroni',
            pts=True,
            n_genes=n_top_genes * 2  # 计算更多基因以便后续筛选
        )

        # 3. 获取所有细胞类型
        cell_types = adata.obs[groupby].unique().tolist()
        cell_types.sort()
        logger.info(f"发现 {len(cell_types)} 个细胞类型: {cell_types}")

        # 4. 为每个细胞类型提取差异基因
        degs_data = {}

        for cell_type in cell_types:
            # 获取该细胞类型的差异表达结果
            result = sc.get.rank_genes_groups_df(adata, group=cell_type)

            # 过滤显著差异基因 (p_adj < 0.05 且 |log2FC| > 0.25)
            significant = result[
                (result['pvals_adj'] < 0.05) &
                (abs(result['logfoldchanges']) > 0.25)
            ]

            if len(significant) == 0:
                logger.warning(f"细胞类型 '{cell_type}' 未找到显著差异基因，使用所有基因")
                significant = result

            # 按log2FC降序排序
            significant = significant.sort_values('logfoldchanges', ascending=False)

            # 取前n_top_genes个基因
            top_genes = significant['names'].head(n_top_genes).tolist()

            # 去除空值和重复
            top_genes = [str(g) for g in top_genes if pd.notna(g)]
            top_genes = list(dict.fromkeys(top_genes))  # 保持顺序去重

            degs_data[cell_type] = top_genes
            logger.info(f"  {cell_type}: {len(top_genes)} 个差异基因")

        # 5. 生成DEGs文件
        # CellPhoneDB官方格式说明:
        #   - 第1列: 基因名称 (每行一个基因)
        #   - 第2列+: 每个细胞类型对应的该基因
        #   - 如果该基因是某细胞类型的DEG，则填入该基因名；否则留空或填"-"

        # 收集所有出现的基因
        all_genes = set()
        for genes in degs_data.values():
            all_genes.update(genes)
        all_genes = sorted(all_genes)

        logger.info(f"总共有 {len(all_genes)} 个差异基因")

        # 创建DataFrame
        # 列1: 基因名称
        # 列2-n: 每个细胞类型中该基因（如果是DEG则填基因名，否则留空）
        degs_columns = ['gene'] + cell_types
        degs_rows = []

        for gene in all_genes:
            row = [gene]
            for cell_type in cell_types:
                # 检查该基因是否是此细胞类型的DEG
                if gene in degs_data[cell_type]:
                    row.append(gene)  # 是DEG，填入基因名
                else:
                    row.append('')    # 不是DEG，留空
            degs_rows.append(row)

        degs_df = pd.DataFrame(degs_rows, columns=degs_columns)

        # 保存文件
        degs_df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"DEGs文件已保存至: {output_path}")
        logger.info(f"DEGs文件形状: {degs_df.shape[0]}行 × {degs_df.shape[1]}列")

        # 打印前几行作为示例
        logger.info("DEGs文件前5行示例:")
        logger.info(f"列名: {list(degs_df.columns)}")
        for i in range(min(5, len(degs_df))):
            row_dict = degs_df.iloc[i].to_dict()
            logger.info(f"  行{i}: gene={row_dict['gene']}, 各细胞类型基因数={[(ct, len(row_dict[ct].split(','))) for ct in cell_types[:3]]}")

        return degs_df

    except Exception as e:
        logger.error(f"生成DEGs文件时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        raise


def export_cellphonedb_input_files(adata, output_dir, groupby='cell_type'):
    """从h5ad文件导出CellPhoneDB所需的计数矩阵和元数据文件"""
    try:
        # 创建输出目录
        cpdb_input_dir = os.path.join(output_dir, "cellphonedb_input")
        os.makedirs(cpdb_input_dir, exist_ok=True)

        # 1. 导出计数矩阵
        counts_path = os.path.join(cpdb_input_dir, "counts.txt")

        # 处理稀疏矩阵 - 大规模数据优化：避免内存爆炸
        if scipy.sparse.issparse(adata.X):
            logger.info("检测到稀疏矩阵，使用稀疏矩阵优化导出...")

            # 分块处理以避免内存溢出
            chunk_size = 10000  # 每次处理10000个细胞
            n_cells = adata.shape[0]
            n_genes = adata.shape[1]

            # 初始化输出文件（写入表头）
            with open(counts_path, 'w') as f:
                # 写入表头：基因名 + 细胞名
                header = '\t'.join([''] + list(adata.obs_names))
                f.write(header + '\n')

            # 分块写入数据
            for chunk_start in range(0, n_cells, chunk_size):
                chunk_end = min(chunk_start + chunk_size, n_cells)
                logger.info(f"处理细胞 {chunk_start+1}-{chunk_end}/{n_cells}...")

                # 提取当前块的稀疏矩阵
                chunk_X = adata.X[chunk_start:chunk_end, :]

                # 转换为密集矩阵（只转换当前块）
                chunk_dense = chunk_X.toarray()

                # 转置为基因 × 细胞格式
                chunk_dense_T = chunk_dense.T

                # 写入当前块
                with open(counts_path, 'a') as f:
                    for gene_idx in range(n_genes):
                        gene_name = adata.var_names[gene_idx]
                        values = '\t'.join(map(str, chunk_dense_T[gene_idx, :]))
                        if chunk_start == 0:
                            f.write(f'{gene_name}\t{values}\n')
                        else:
                            f.write(f'{values}\n')

            logger.info(f"计数矩阵已导出至: {counts_path}")
        else:
            # 非稀疏矩阵，直接处理
            counts_df = pd.DataFrame(adata.X.T, index=adata.var_names, columns=adata.obs_names)
            counts_df.to_csv(counts_path, sep='\t')
            logger.info(f"计数矩阵已导出至: {counts_path}")

        # 2. 导出元数据
        meta_path = os.path.join(cpdb_input_dir, "meta.txt")
        meta_df = adata.obs[[groupby]].copy()
        meta_df.index.name = 'Cell'
        meta_df.columns = ['cell_type']
        meta_df.to_csv(meta_path, sep='\t')
        logger.info(f"元数据已导出至: {meta_path}")

        return counts_path, meta_path

    except Exception as e:
        logger.error(f"导出CellPhoneDB输入文件时出错: {str(e)}")
        raise


def validate_degs_file(degs_path):
    """验证DEGs文件格式是否符合CellPhoneDB标准"""
    try:
        logger.info(f"验证DEGs文件: {degs_path}")

        degs_df = pd.read_csv(degs_path, sep='\t')

        # 检查必需列
        if len(degs_df.columns) < 2:
            logger.error(f"DEGs文件格式错误: 列数不足 ({len(degs_df.columns)}列)")
            return False

        # 检查第一列是否为基因列
        first_col = degs_df.columns[0]
        if first_col not in ['gene', 'Gene', 'gene_name', 'Gene_name']:
            logger.warning(f"DEGs文件第一列名称为 '{first_col}'，可能不是基因列")

        # 检查数据内容
        if degs_df.shape[0] < 10:
            logger.error(f"DEGs文件基因数量过少: {degs_df.shape[0]}")
            return False

        logger.info(f"DEGs文件验证通过: {degs_df.shape[0]} 个基因, {degs_df.shape[1]-1} 个细胞类型")
        return True

    except Exception as e:
        logger.error(f"验证DEGs文件时出错: {str(e)}")
        return False


def run_cellphonedb_analysis(h5ad_path, cpdb_zip, degs_path, microenvs_path,
                             groupby='cell_type', output_dir=None, generate_degs=False,
                             degs_method='wilcoxon', n_top_genes=250):
    """
    使用DEGs和微环境文件运行CellPhoneDB分析

    参数:
        h5ad_path: h5ad文件路径
        cpdb_zip: CellPhoneDB数据库zip文件路径
        degs_path: 差异表达基因文件路径
        microenvs_path: 微环境文件路径
        groupby: 分组列名
        output_dir: 输出目录
        generate_degs: 是否自动生成DEGs文件
        degs_method: DEGs计算方法
        n_top_genes: 每个细胞类型的差异基因数量
    """
    try:
        # 设置输出目录
        if output_dir is None:
            output_dir = os.path.dirname(h5ad_path)

        # 1. 读取h5ad文件
        logger.info(f"读取h5ad文件: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)

        # 检查分组列是否存在
        if groupby not in adata.obs.columns:
            raise ValueError(f"分组列 '{groupby}' 不存在于adata.obs中")

        logger.info(f"数据集信息: {adata.shape[0]}个细胞, {adata.shape[1]}个基因")

        # 2. (可选) 生成DEGs文件
        if generate_degs:
            logger.info("生成DEGs文件...")
            degs_path = os.path.join(output_dir, "DEGs.tsv")
            generate_degs_file(adata, degs_path, groupby, degs_method, n_top_genes)
        else:
            # 验证DEGs文件格式
            if not validate_degs_file(degs_path):
                logger.error("DEGs文件验证失败")
                return False

        # 3. 导出CellPhoneDB所需的计数矩阵和元数据
        logger.info("导出CellPhoneDB输入文件...")
        counts_path, meta_path = export_cellphonedb_input_files(
            adata, output_dir, groupby
        )

        # 4. 检查所有必需文件是否存在
        logger.info("检查所有必需文件是否存在...")
        required_files = [cpdb_zip, counts_path, meta_path, degs_path, microenvs_path]
        for file in required_files:
            if not os.path.exists(file):
                logger.error(f"错误: 文件不存在 - {file}")
                return False

        # 5. 验证微环境文件格式
        try:
            logger.info("验证微环境文件格式...")
            microenvs_df = pd.read_csv(microenvs_path, sep='\t', header=None)
            if len(microenvs_df.columns) < 2:
                logger.error("微环境文件格式错误：列数不足")
                return False
            logger.info("微环境文件格式验证通过")
        except Exception as e:
            logger.error(f"验证微环境文件时出错: {str(e)}")
            return False

        # 7. 设置CellPhoneDB结果目录（使用动态时间戳）
        from datetime import datetime
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_dir = os.path.join(output_dir, f"cellphonedb_results_{timestamp}")
        os.makedirs(results_dir, exist_ok=True)

        # 8. 智能检测基因标识符类型
        try:
            counts_df_check = pd.read_csv(counts_path, sep='\t', nrows=5)
            first_gene_name = str(counts_df_check.iloc[0, 0])

            if first_gene_name.startswith('ENSG'):
                possible_count_types = ['ensembl', 'gene_name', 'hgnc_symbol']
                logger.info(f"检测到 Ensembl 基因名格式")
            elif '-' in first_gene_name and len(first_gene_name.split('-')) > 1:
                possible_count_types = ['gene_name', 'hgnc_symbol', 'ensembl']
                logger.info(f"检测到复杂基因名格式")
            else:
                possible_count_types = ['gene_name', 'hgnc_symbol', 'ensembl']
                logger.info(f"检测到标准基因名格式")

            logger.info(f"示例基因名: {first_gene_name}")
        except Exception as e:
            logger.warning(f"无法检测基因名格式，使用默认顺序: {str(e)}")
            possible_count_types = ['gene_name', 'hgnc_symbol', 'ensembl']

        # 9. 运行CellPhoneDB分析（尝试多种基因标识符）
        success = False
        for count_type in possible_count_types:
            try:
                logger.info(f"尝试使用计数数据类型: {count_type}")
                logger.info("开始CellPhoneDB分析...")

                cpdb_degs_analysis_method.call(
                    cpdb_file_path=cpdb_zip,
                    meta_file_path=meta_path,
                    counts_file_path=counts_path,
                    degs_file_path=degs_path,
                    counts_data=count_type,
                    output_path=results_dir,
                    subsampling=False,           # 大规模数据不使用子采样
                    iterations=1000,             # 迭代次数
                    threshold=0.05,              # 显著性阈值
                    threads=multiprocessing.cpu_count(),
                    result_precision=3
                )

                logger.info("✓ CellPhoneDB分析成功完成!")
                success = True
                break

            except Exception as e:
                logger.warning(f"✗ 使用 '{count_type}' 时分析失败: {str(e)}")
                logger.warning("尝试下一种基因标识符类型...")

        if not success:
            logger.error("所有基因标识符类型尝试均失败")
            return False

        logger.info(f"分析完成! 结果保存在: {results_dir}")
        logger.info("重要结果文件:")
        logger.info(f"  - deconvoluted.txt: 细胞对互作详细结果")
        logger.info(f"  - means.txt: 细胞类型对平均互作分数")
        logger.info(f"  - pvalues.txt: 统计显著性p值")
        logger.info(f"  - significant_means.txt: 显著互作汇总")

        return True

    except Exception as e:
        logger.error(f"运行CellPhoneDB分析时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False


def prepare_degs_only(h5ad_path, output_path, groupby='cell_type',
                      method='wilcoxon', n_top_genes=250):
    """
    仅生成DEGs文件，不运行CellPhoneDB分析

    参数:
        h5ad_path: h5ad文件路径
        output_path: 输出DEGs文件路径
        groupby: 分组列名
        method: 差异分析方法
        n_top_genes: 每个细胞类型的差异基因数量
    """
    try:
        logger.info("读取h5ad文件...")
        adata = sc.read_h5ad(h5ad_path)

        if groupby not in adata.obs.columns:
            raise ValueError(f"分组列 '{groupby}' 不存在于adata.obs中")

        logger.info(f"数据集信息: {adata.shape[0]}个细胞, {adata.shape[1]}个基因")

        # 生成DEGs文件
        degs_df = generate_degs_file(adata, output_path, groupby, method, n_top_genes)

        logger.info("DEGs文件生成完成!")
        return True

    except Exception as e:
        logger.error(f"生成DEGs文件时出错: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
        return False


if __name__ == '__main__':
    multiprocessing.freeze_support()

    # 设置命令行参数
    parser = argparse.ArgumentParser(
        description='从h5ad数据生成DEGs文件并运行CellPhoneDB分析'
    )
    parser.add_argument('--mode', choices=['full', 'degs_only'], default='full',
                        help='运行模式: full=完整流程, degs_only=仅生成DEGs')
    parser.add_argument('--h5ad', type=str,
                        default=r"G:\DATA\SingleCell\TabulaSapiens\tabula_sapiens_liver\c264e09f-7c3b-4294-b0f4-82a790bd0014.h5ad",
                        help='输入h5ad文件路径')
    parser.add_argument('--cpdb', type=str, default="./v5.0.0/cellphonedb.zip",
                        help='CellPhoneDB数据库zip文件路径')
    parser.add_argument('--degs', type=str, default=r"./liver_DEGs_wilcoxon.tsv",
                        help='差异表达基因文件路径 (full模式下如果--generate_degs为True则忽略)')
    parser.add_argument('--microenv', type=str, default="liver_microenvironment.tsv",
                        help='微环境文件路径')
    parser.add_argument('--groupby', type=str, default='cell_type',
                        help='分组列名 (默认: cell_type)')
    parser.add_argument('--outdir', type=str, default="./",
                        help='输出目录 (默认: h5ad文件所在目录)')
    parser.add_argument('--generate_degs', action='store_true',
                        help='是否自动生成DEGs文件 (full模式下)')
    parser.add_argument('--degs_method', type=str, default='wilcoxon',
                        choices=['wilcoxon', 't-test', 'logreg'],
                        help='DEGs计算方法 (默认: wilcoxon)')
    parser.add_argument('--n_top_genes', type=int, default=250,
                        help='每个细胞类型的差异基因数量 (默认: 250)')

    args = parser.parse_args()

    logger.info("=== 开始细胞互作分析流程 ===")
    logger.info(f"运行模式: {args.mode}")

    try:
        if args.mode == 'degs_only':
            # 仅生成DEGs文件
            degs_output = args.outdir
            if os.path.isdir(degs_output):
                degs_output = os.path.join(degs_output, "DEGs.tsv")

            success = prepare_degs_only(
                h5ad_path=args.h5ad,
                output_path=degs_output,
                groupby=args.groupby,
                method=args.degs_method,
                n_top_genes=args.n_top_genes
            )

            if success:
                logger.info("DEGs文件生成完成!")
            else:
                logger.error("DEGs文件生成失败")

        else:  # full模式
            success = run_cellphonedb_analysis(
                h5ad_path=args.h5ad,
                cpdb_zip=args.cpdb,
                degs_path=args.degs,
                microenvs_path=args.microenv,
                groupby=args.groupby,
                output_dir=args.outdir,
                generate_degs=args.generate_degs,
                degs_method=args.degs_method,
                n_top_genes=args.n_top_genes
            )

            if success:
                logger.info("细胞互作分析成功完成!")
            else:
                logger.error("细胞互作分析失败")

    except Exception as e:
        logger.error(f"分析过程中发生错误: {str(e)}")
        import traceback
        logger.error(traceback.format_exc())
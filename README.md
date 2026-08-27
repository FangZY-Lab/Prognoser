# Prognoser

`Prognoser` 是一个整合多组表达谱与生存信息、发现预后基因标签（prognostic gene
signatures）的 R 包。它在每个数据集中进行单变量 Cox 回归，综合多个数据集筛选保护型
（HR < 1）和风险型（HR >= 1）基因，并通过聚类去除冗余，最终输出非冗余基因集。

## 依赖

运行前请确认已安装以下 R 包：

- `survival`
- `RobustRankAggreg`
- `meta`

```r
install.packages(c("survival", "RobustRankAggreg", "meta"))
```

## 安装

### 方式一：从源码压缩包安装

```r
install.packages("~/Desktop/Prognoser_0.1.0.tar.gz", repos = NULL, type = "source")
```

### 方式二：从 GitHub 安装

```r
install.packages("remotes")
remotes::install_github("FangZY-Lab/Prognoser")
```

安装完成后加载：

```r
library(Prognoser)
```

## 数据格式

`Prognoser()` 通过 `get()` 从全局环境中读取数据，要求：

- 每个数据集是一个 `data.frame`，行名为基因名，列名为样本名。
- 每个数据集对应一个生存数据框，命名规则为 `<数据集名>_survival`，包含 `time`
  （数值）和 `status`（数值，`0` = 存活，`1` = 死亡）两列，行名为与表达矩阵列名
  一致的样本名。

## 生成模拟数据

下面的代码生成两个模拟数据集（`SimDataA`、`SimDataB`）及其生存数据。为了演示，
在表达矩阵中植入了 10 个“风险基因”（`GENE1`–`GENE10`，表达越高生存越差）和
10 个“保护基因”（`GENE11`–`GENE20`，表达越高生存越好）。

```r
set.seed(123)

make_dataset <- function(n_samples = 80, n_genes = 100, seed = 1) {
  set.seed(seed)

  # 潜在风险分数：分数越高，生存越差
  risk_score <- rnorm(n_samples)

  # 生存数据：风险越高，生存时间越短
  time <- rweibull(n_samples, shape = 1, scale = exp(2 - risk_score))
  status <- rbinom(n_samples, 1, 0.7)
  surv <- data.frame(
    time = round(time, 3),
    status = status,
    row.names = paste0("S", seq_len(n_samples))
  )

  # 表达矩阵：行为基因，列为样本
  mat <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)

  # 前 10 个基因为风险基因（与 risk_score 正相关）
  mat[1:10, ] <- mat[1:10, ] + matrix(rep(risk_score * 2, each = 10), nrow = 10)

  # 第 11-20 个基因为保护基因（与 risk_score 负相关）
  mat[11:20, ] <- mat[11:20, ] - matrix(rep(risk_score * 2, each = 10), nrow = 10)

  expr <- as.data.frame(mat)
  rownames(expr) <- paste0("GENE", seq_len(n_genes))
  colnames(expr) <- paste0("S", seq_len(n_samples))

  list(expr = expr, surv = surv)
}

dA <- make_dataset(seed = 11)
dB <- make_dataset(seed = 22)

# 放入全局环境（函数通过 get() 读取）
SimDataA <- dA$expr
SimDataA_survival <- dA$surv
SimDataB <- dB$expr
SimDataB_survival <- dB$surv
```

## 运行

```r
# 先验 / 知识基因集（可选，作为候选基因池的一部分）
knowledge <- list(immune = c("GENE1", "GENE2", "GENE3", "GENE11", "GENE12"))

result <- Prognoser(
  expression_accession_vector = c("SimDataA", "SimDataB"),
  knowledge_genesets = knowledge,
  denovo_threshold = 0.05,
  gene_threshold = 0.01,
  min_genes = 3,
  similarity_threshold = 0.5,
  linkage = "ward.D"
)

# 最终的非冗余基因集
str(result$genesets)

# 保护方向 / 风险方向的综合效应矩阵
head(result$combined_matrix_protection)
head(result$combined_matrix_risk)
```

示例中使用 `gene_threshold = 0.01`、`min_genes = 3` 是为了让演示结果更稳定；实际
分析时可按需求调整，默认值见下表。

## 参数说明

| 参数 | 说明 | 默认值 |
| --- | --- | --- |
| `expression_accession_vector` | 数据集名向量，每个数据集需在全局环境中存在，并配套 `<名>_survival` | 无 |
| `denovo_threshold` | 单变量 Cox 的 P 值阈值（筛选保护 / 风险基因） | `0.05` |
| `knowledge_genesets` | 先验基因集列表（由字符向量组成的 `list`） | 无 |
| `na_ratio` | 基因在数据集间缺失比例的上限 | `0.3` |
| `gene_threshold` | 稳健秩聚合（RRA）的 P 值阈值 | `0.001` |
| `min_genes` | 最终基因集的最小基因数 | `5` |
| `similarity_threshold` | 基因集聚类的 Jaccard 相似度阈值 | `0.5` |
| `linkage` | 层次聚类方法，需传入单个值（如 `"ward.D"`） | 见下方注意 |

## 返回结果

- `genesets`：最终的非冗余基因集列表。低风险（保护）基因集以 `LRGS` 前缀命名，
  高风险基因集以 `HRGS` 前缀命名。
- `combined_matrix_protection`：保护方向基因的综合效应矩阵（符号 × `-log10(P)`），
  含 `p_protection` 列。
- `combined_matrix_risk`：风险方向基因的综合效应矩阵，含 `p_risk` 列。

## 注意

1. **`linkage` 参数必须显式传入单个方法**（如 `"ward.D"`、`"average"`、
   `"complete"`）。如果不显式指定，会使用原函数的向量默认值，导致 `hclust()`
   报错。
2. 数据必须放在全局环境，且命名严格遵循 `<数据集名>` 与 `<数据集名>_survival`。

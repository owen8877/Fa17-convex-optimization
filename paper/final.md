# 凸优化大作业报告

__陈子恒 1500010632__

## 第一题

我们已经非常熟悉有关于离散最佳传输问题（Discrete Optimal Transport，下简称 __DOT__ ）的线性规划表示了：
$$
\begin{align}
\textbf{min} & \sum_{i=1}^m \sum_{j=1}^n c_{ij} \pi_{ij} \\
\textbf{s.t.} & \sum_{j=1}^n \pi_{ij}=\mu_i, \forall i=1, \dots m \\ 
& \sum_{i=1}^m \pi_{ij}=\nu_j, \forall j=1, \dots n \\ 
& \pi_{ij} \ge 0.
\end{align}
$$
注意到这里传输方案 $\pi$ 是矩阵形式，写成向量形式就是
$$
\begin{align}
\textbf{min} \quad & f^T x\\
\textbf{s.t.} \quad & A x=b \\
&x \ge 0
\end{align}
$$
这里
$$
\begin{align}
f &= (c_{11}, c_{21}, \dots c_{n1}, \dots, c_{1n}, \dots c_{nn})^T \\
x &= (\pi_{11}, \pi_{21}, \dots \pi_{n1}, \dots, \pi_{1n}, \dots \pi_{nn})^T \\
A &= \begin{pmatrix} I_m & \cdots & I_m \\ E_1 &\cdots & E_n \end{pmatrix} , (E_j)_{st}=\delta_{sj} \\
b &= (\mu_1, \dots \mu_m, \nu_1, \dots \nu_n)^T.
\end{align}
$$
根据此，我们可以使用 mosek 的 `linprog` 函数计算 DOT （详见 `simple-tests/dir-mosek.m`）：

```matlab
switch opts.method
    case 'interior'
        options.Interior = 'on';
    case 'simplex'
        options.Simplex = 'on';
    otherwise
        options.Interior = 'on';
end
options.Diagnostics = 'on';
options.Display = 'iter';

[x, ~, ~, ~, ~] = linprog(f, A, b, B, c, l, u, x0, options);
```

我们分别对内点法和单纯性法进行测试，这里 $m=n=64$ ：

| method             | cpu     | error    | funcval  |
| ------------------ | ------- | -------- | -------- |
| dir_mosek_interior | 0.46802 | 0        | 0.023837 |
| dir_mosek_simplex  | 0.28712 | 5.67e-12 | 0.023837 |

这表明单纯性法较内点法有较大的优势，这也是下文所着重要讨论的问题。



## 第二题

既然 DOT 是一个有等式约束的优化问题，我们可以对原问题做 ADMM ，推导如下：

> 引入变量 $z$ ，将不等式约束隐含在目标函数中：
> $$
> \begin{align}
> \textbf{min} \quad & f^T x + I_C(z)\\
> \textbf{s.t.} \quad & \begin{pmatrix} A \\ I\end{pmatrix} x + \begin{pmatrix} 0 \\ -I\end{pmatrix}z=\begin{pmatrix} b \\ 0\end{pmatrix}, C=\{z_i \ge 0, \forall i\}
> \end{align}
> $$
> 引入乘子 $y_1 \in \mathbf{R}^{m+n}$ 与 $y_2 \in \mathbf{R}^{mn}$ ，得到迭代格式
> $$
> \begin{align}
> v_1 &\leftarrow b+y^{(k)}_1/\beta \\
> v_2 &\leftarrow z+y^{(k)}_1/\beta \\
> x^{(k+1)} &\leftarrow (A^TA+I) \backslash (A^Tv_1+v_2-f/\beta) \\
> z^{(k+1)} &\leftarrow P_C(x^{(k+1)}-y^{(k)}_2/\beta) \\
> y_1^{(k+1)} &\leftarrow y_1^{(k)} + (b-Ax^{(k+1)}) \\
> y_2^{(k+1)} &\leftarrow y_2^{(k)} + (z^{(k+1)}-x^{(k+1)})
> \end{align}
> $$
> 第三步左除是用 CG 完成的，否则精度很差。

ADMM 效果并不佳；以下是 $m=n=64$ 的测试结果：

| method             | cpu     | error      | funcval  |
| ------------------ | ------- | ---------- | -------- |
| dir_mosek_interior | 0.38076 | 0          | 0.028047 |
| dir_mosek_simplex  | 0.35575 | 1.9222e-08 | 0.028047 |
| admm_nesterov      | 8.5171  | 0.00045024 | 0.028052 |
| admm               | 6.9257  | 0.0023568  | 0.027985 |

一方面误差很大，这是由于等式约束难以起到效果的缘故。另一方面是所需时间远长于 mosek 的通用算法，这是由于 ADMM 算法会在 $x_i \approx \frac{1}{mn}$ 处停留很长时间，直至某几个变量开始出现明显的增长才使收敛步骤开始；当然这也与这是一阶算法有关。

对于对偶问题当然也可以做 ADMM ，推导步骤如下：

> 对偶问题是
> $$
> \begin{align}
> \textbf{min} \quad & b^T \nu\\
> \textbf{s.t.} \quad & f+A^T\nu \ge 0
> \end{align}
> $$
> 改写成 ADMM 的形式
> $$
> \begin{align}
> \textbf{min} \quad & b^T \nu + I_C(\lambda)\\
> \textbf{s.t.} \quad & A^T\nu - \lambda = -f
> \end{align}
> $$
> 于是迭代格式为
> $$
> \begin{align}
> \nu^{(k+1)} &\leftarrow (A^TA+\epsilon I) \backslash \left(\frac{Ax^{(k)}-b}{\beta}-A(f-\lambda^{(k)})\right) \\
> \lambda^{(k+1)} &\leftarrow P_C(A^T \nu^{(k+1)} +f - x^{(k)} /\beta) \\
> x^{(k+1)} &\leftarrow x^{(k)} - (f + A^T \nu^{(k+1)}-\lambda^{(k+1)})
> \end{align}
> $$
> 更新 $\nu$ 的一步中 $\epsilon$ 的作用是为了使左除数值稳定。

算法已经在 `admm_dual.m` 中实现了，但是难以在合理的时间内收敛。



## 第三题

### Transport Simplex Method

本质上来说以下所实现的三个算法都是组合优化算法，这是因为单纯形法——进一步的，传输单纯性法，充分利用了 $A$ 的稀疏结构。由于它们充分考虑了 $A$ 的特殊性，速度是很快的（以下测试为 $m=n=64$ ）：

| method                 | cpu     | error      | funcval  |
| ---------------------- | ------- | ---------- | -------- |
| dir_mosek_interior     | 0.38611 | 0          | 0.027697 |
| dir_mosek_simplex      | 0.36667 | 4.5602e-10 | 0.027697 |
| admm_nesterov          | 10.703  | 0.00083133 | 0.027701 |
| admm                   | 7.219   | 0.0025004  | 0.027714 |
| transimplex            | 0.27942 | 4.5602e-10 | 0.027697 |
| tran_simplex_normal    | 0.56725 | 4.5602e-10 | 0.027697 |
| tran_simplex_improved  | 0.64472 | 4.5602e-10 | 0.027697 |
| tran_simplex_shortlist | 0.78707 | 4.5602e-10 | 0.027697 |

完全使用 matlab 实现的单纯性法 transimplex 速度就已经能与 mosek 相提并论。

但是当规模增长后开始迅速变慢（以下测试为 $m=n=256$）：

| method                 | cpu     | error      | funcval   |
| ---------------------- | ------- | ---------- | --------- |
| dir_mosek_interior     | 0.73828 | 0          | 0.0067589 |
| dir_mosek_simplex      | 0.50615 | 2.9584e-09 | 0.0067589 |
| transimplex            | 4.6906  | 2.9584e-09 | 0.0067589 |
| tran_simplex_normal    | 0.59002 | 2.9584e-09 | 0.0067589 |
| tran_simplex_improved  | 0.65727 | 2.9584e-09 | 0.0067589 |
| tran_simplex_shortlist | 0.75122 | 2.9584e-09 | 0.0067589 |

我们使用 matlab 的 profile 工具得知主要时间花在了利用深度优先搜索确定被剔除变量的位置上，所以我们改用 C++ 和 MEX API 实现算法；结果表明速度是很快的。

#### 算法分析

传输单纯形法求解主要有两个步骤：构造初始解和更新解。

##### 初始解的构造

根据单纯性法，由于我们有 $m+n-1$ 个约束条件，所以其实只有 $m+n-1$ 个变量是自有变量，剩下的都是约束变量。根据单纯性法，如果选择 $x_{k_1}, \dots x_{k_{m+n-1}}$ 作为自由变量能使得目标函数是它们的非正线性组合，我们就达到了最优解。

由于只有自由变量能够取零，因此一个自然的想法是我们挑选每行/列代价最大的变量作为自由变量，并且尽可能满足这个变量的传输直至这一行或列的总和用尽；这就是最小费用法则（Least Cost Rule）：
$$
\begin{align}
& searchOnRow \leftarrow true \\
& row \leftarrow 1 \\
& column \leftarrow 1 \\
& \textbf{for} \quad i=1:(m+n-1) \\
& \qquad \textbf{if} \quad searchOnRow \\
& \qquad \qquad column \leftarrow \textbf{argmax}_l \quad c_{row, l} \\
& \qquad \qquad \textbf{if} \quad \mu_{row, column} \ge \nu_{row, column} \\
& \qquad \qquad \qquad \mu_{row, column} \leftarrow \mu_{row, column} - \nu_{row, column} \\
& \qquad \qquad \textbf{else} \\
& \qquad \qquad \qquad \nu_{row, column} \leftarrow \nu_{row, column} - \mu_{row, column} \\
& \qquad \qquad \qquad searchOnRow \leftarrow false \\
& \qquad \textbf{else} \\
& \qquad \qquad row \leftarrow \textbf{argmax}_l \quad c_{l, column} \\
& \qquad \qquad \textbf{if} \quad \nu_{row, column} \ge \mu_{row, column} \\
& \qquad \qquad \qquad \nu_{row, column} \leftarrow \nu_{row, column} - \mu_{row, column} \\
& \qquad \qquad \textbf{else} \\
& \qquad \qquad \qquad \mu_{row, column} \leftarrow \mu_{row, column} - \nu_{row, column} \\
& \qquad \qquad \qquad searchOnRow \leftarrow true \\
& \qquad \rm{put \quad (row, column) \quad into \quad query} \\
& \qquad c_{row, column} \leftarrow \infty \\
\end{align}
$$

##### 解的优化

为了判断一个解是不是最优的，我们先要构造对偶变量 $u \in \mathbf{R}^{m\times 1},v\in \mathbf{R}^{1\times n}$ ，并根据相对价格 $c - u \mathbf{1}^T - \mathbf{1}v^T$ 是否有负元素来判断一个解是否达到最优。详见 10.1007 。

在找出有负相对价格的变量后，可以断定这个变量的加入使得自由变量间以同行/列为边的图产生了有且仅有一个环，那么，找到这个环上距离新加入变量距离为偶数步且运输量最少的元素，就是这一次更新需要更新的运输量。详见 `normaltransimplex.cpp` 。

####  改进要点

除了算法步骤必要的优化外，在每一步迭代的时候可以存储许多关键的信息避免重复构造

- 节点的相互连接使用图存储
- 初始解的选取（即最小费用原则）
- 在判断解是否是最优的时候，我们不需要遍历整个列表寻找最小相对价格元素；经验估计一般直接更新第一个就可以了，虽然迭代步数较多，但是搜索时间显著变短。
- 接上一条，在遍历非自由变量表时可以从上一次搜索的位置的下一行继续搜索，这样能够在较短时间内更新较多变量。见 `minimalrowtransimplex.cpp` 。

### Shortlist Method

在 pone 011 中作者详细比较了不同初始解算法是否能给出一个较好的解。 Shortlist 算法的思想是预先储存一些和这个元素费用较小的邻居，在更新的时候优先搜索这些邻居；按理说这个思想是很好的，但是实际测试的时候效果不显著，主要是因为如果是随机数据的话 shortlist 基本没有价值，而如果是欧式距离的话速度肯定也没有 shielding 快。

我们开始使用图片集作为输入数据。


| method                 | cpu     | error   | funcval | average cpu |
| ---------------------- | ------- | ------- | ------- | ----------- |
| dir_mosek_simplex      | 2.3936  | 0       | 610.69  | 2.3936      |
| tran_simplex_normal    | 1.8187  | 0.46795 | 610.69  | 1.8187      |
| tran_simplex_improved  | 1.6048  | 0.49995 | 610.69  | 1.6048      |
| tran_simplex_shielding | 0.87166 | 0.48548 | 610.69  | 0.87166     |
| tran_simplex_shortlist | 1.8504  | 0.53429 | 610.69  | 1.8504      |

#### 图片的生成

根据 0778 ，我们一共需要生成 10 类数据，目前我们仅生成了第一种数据，分辨率从 $8\times 8$ 到 $48\times 48$ 不等。上一节的数据是 $ 16 \times 16$ 分辨率的白噪声图片在 $p=2$ 的欧式距离代价。

### Shielding Method

这个算法的思想是如果我们能够给每一个元素 $x$ 找一族“邻居” $\Lambda_x$ ，使得
$$
\exists x^* \in \Lambda_x \quad \textbf{s.t.} \quad Cost(\tilde{x}, x) > Cost(\tilde{x}, x^*) + Cost(x^*, x)
$$
那么我们只需要每次检验邻居的相对费用是不是负的就可以了。在最简单的 $p=2$ 情形（即 $Cost(x, y) = ||x-y||^p$）下， $\Lambda_x$ 就是 $x$ 周围的八个邻居。（有没有理论保证收敛速度来着？）

### AHA Method

所限于时间问题我们无法给出这个算法的实现。
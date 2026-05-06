<h1 align = "center">矩阵理论与应用II</h1>
<h3 align="right">习题求解整理</h3>
<p align="right">SY2502110-朱宛瑜<br>Date: 2026.5.6</p>

# 1.自测题

1.  设 \(A, B \in \mathbb{C}^{n \times n}\)，且 \(A\) 的特征值全为实数，证明：矩阵方程 \(X + AXA + A^2 XA^2 = B\) 有唯一解。

2.  证明：\(A \in \mathbb{R}_r^{m \times n}\)，\(\sigma_i\) \((1 \leq i \leq r)\) 为奇异值，则 \(\operatorname{tr}(A^{\mathrm{H}}A) = \sum_{i = 1}^{r} \sigma_i^2\)。

3.  求下列矩阵的奇异值分解：

    

4.  设 \(A\) 是正规矩阵，证明：
    (1) \(A^{+}A = AA^{+}\)；
    (2) \((A^{+})^{k} = (A^{k})^{+}\)（\(k\) 为正整数）。

5.  验证方程是相容的，并用 \(A^{(1)}\) 表示通解。

    

6.  求上题通解中的极小范数解。

7.  求解矩阵方程 \(AX + XB = C\)，其中

    

8.  设

\[
A = \begin{pmatrix}
1 & 0 & -1 & 1 \\
0 & 2 & 2 & 2 \\
-1 & 4 & 5 & 3
\end{pmatrix},\quad
b = \begin{pmatrix}
4 \\ -2 \\ -2
\end{pmatrix},
\]
(1) 用广义逆矩阵方法判定方程组 \(Ax = b\) 是否相容；

(2) 求 \(Ax = b\) 的极小范数解或极小最小二乘解。

9.  设
    \[
    A(t) = \begin{pmatrix}
    e^{2t} & t e^{t} \\
    0 & e^{-t}
    \end{pmatrix},
    \]
    (1) 证明 \(A(t)\) 可逆并求 \(A^{-1}(t)\)；
    (2) 求 \(\dfrac{\mathrm{d}}{\mathrm{d}t} A^{-1}(t)\)。

10. 设
   \[
   A(t) = \begin{pmatrix}
   t \cos t & e^{t} \sin t \\
   t^{2} + 1 & \ln (1 + t)
   \end{pmatrix},
   \]
   计算：
   (1) \(\dfrac{\mathrm{d}A(t)}{\mathrm{d}t}\)；
   (2) \(\displaystyle\int_{0}^{t} A(s) \mathrm{d}s\)。

11. 设
    \[
    A = \begin{pmatrix}
    -1 & -2 & 6 \\
    -1 & 0 & 3 \\
    -1 & -1 & 4
    \end{pmatrix},\quad
    b(t) = \begin{pmatrix}
    e^{t} \\ 2e^{t} \\ e^{t}
    \end{pmatrix},
    \]
    求 \(e^{At}\)，并求微分方程
    \[
    \frac{\mathrm{d}x(t)}{\mathrm{d}t} = A x(t) + b(t)
    \]
    满足初始条件 \(x(0) = (-1,0,1)^{\mathrm{T}}\) 的解。

12. 证明：若 \(AX = b\) 不相容，则 \(A^{H} A X = A^{H} b\) 一定有解，且它的解也是 \(AX = b\) 的最小二乘解。

# 2.往年题
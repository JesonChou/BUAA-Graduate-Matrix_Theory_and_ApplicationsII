<h1 align = "center">矩阵理论与应用II复习</h1>
<p align="right">SY2502110 朱宛瑜</p>

# 1.基本概念

## 1.1 广义逆矩阵及其应用

### 1.1.1 基本概念及奇异值分解

1. 定义1.1.1（广义逆矩阵）：设$A\in \mathbb{C}^{m×n}$，若$X\in \mathbb{C}^{n×m}$满足以下四个Penrose方程的全部或一部分：

   > [!IMPORTANT]
   >
   > 1、$AXA=A$
   >
   > 2、$XAX=X$
   >
   > 3、$(AX)^H=AX$
   >
   > 4、$(XA)^H=XA$

   则称矩阵$X$是矩阵$A$的==广义逆矩阵==。

   符号表示为$G=A^{(i,j,k)}$，其中$i,j,k=1,2,3,4$且$i,j,k$互不相等。

   满足全部四条的记为$G=A^{(1,2,3,4)}$或者$A^+$，称为加号逆、伪逆或者$Moore-Penrose$广义逆。

   $A$为奇异矩阵或者非方阵时，$A$的减号逆矩阵不一定唯一，记作集合$A\{1\}$。其他逆同样。

2. 常见广义逆为$A^{(1)}$，$A^{(1,3)}$，$A^{(1,4)}$和$A^+$。其中$A^{(1)}$称为减号逆，记作$A^-$；$A^{(1,3)}$称为最小二乘广义逆，记作$A^-_l$；$A^{(1,4)}$称为极小范数广义逆，记为$A_m^-$。

3. $Jordan$标准形分解：$J=P^-AP=diag(J_1,J_2,...,J_s)$。

4. $Jordan$标准形分解的局限性：

   > [!IMPORTANT]
   >
   > 1、矩阵$A$必须是方阵；
   >
   > 2、$Jordan$阵不一定是对角阵；
   >
   > 3、相似变换矩阵不是酉矩阵。
   >
   > 注：酉矩阵的定义：$A^HA=AA^H=I$。即复空间上的正交矩阵（$A^TA=AA^T=I$）。

   但是对于$A\in\mathbb{C}^{m×n}$，$m$不一定等于$n$，即$A$是==非方阵==。能否找到一种可以克服$Jordan$分解的局限性的分解方法？

5. 基于上述思路，我们选择Hermite矩阵（即$A^H=A$）$A^HA_{n×n}$与$AA^H_{m×m}$，则分别存在酉矩阵$U\in\mathbb{C}^{m×m}$与$V\in\mathbb{C}^{n×n}$，使得
   $$
   U^HAA^HU=\Sigma_1,V^HA^HAV=\Sigma_2
   $$
   式中$\Sigma_1$与$\Sigma_2$为对角矩阵，对角线元素为矩阵$AA^H$和$A^HA$的==特征值==。<font color='red'>（这里涉及到Hermite矩阵的性质，Hermite矩阵是正规矩阵，其一定可以被酉相似对角化）</font>

   则矩阵$A$可以分解为
   $$
   A=UDV^H
   $$
   其中$DD^H=\Sigma_1$，$D^HD=\Sigma_2$。

   > [!IMPORTANT]
   >
   > 具体推导步骤为：
   > $$
   > AA^H=UDV^H·(UDV^H)^H=UDV^HVD^HU^H
   > $$
   > 由于$U$和$V$均为酉矩阵，则根据酉矩阵的性质（$A^HA=AA^H=I$）可得
   > $$
   > AA^H=UDD^HU^H
   > $$
   > 将式$(4)$中的矩阵左乘$U^H$并右乘$U$（因为$U$是酉矩阵，所以$U^HU=UU^H=I$）：
   > $$
   > U^HAA^HU=DD^H=\Sigma_1
   > $$
   > 对于$V$和$A^HA$同理。

6. ==引理1.1.1== 设$A\in\mathbb{C}^{m×n}$，则$rank(A^HA)=rank(AA^H)=rank(A)$。

   > [!NOTE]
   >
   > （以下推导过程需要了解掌握，可能出证明题）
   >
   > 首先需要知道矩阵秩的性质，即$rank(A)=rank(A^H)$。
   >
   > 设$\vec x\in N(A^HA)$，则有$A^HA\vec{x}=0$。

7. ==引理1.1.2== 设$A\in\mathbb{C}^{m×n}$，则$A^HA$与$AA^H$都是半正定Hermite矩阵。

   > [!NOTE]
   >
   > 矩阵$A$是半正定Hermite矩阵还有以下命题等价（PPT P13）
   >

8. ==引理1.1.3== 设$A\in\mathbb{C}^{m×n}_r$，则$A^HA$与$AA^H$的所有非零特征值完全相同且**非零特征值的个数**均为$r$。<font color='red'>（这里与奇异值的定义相关，$A^HA$与$AA^H$的非零特征值个数为$r$，则矩阵$A$的正奇异值的个数为$r$）</font>

9. 计算$A^HA$和$AA^H$的特征值时，可优先考虑计算阶数较小的矩阵。

10. ==定义1.1.2（奇异值）== 设$A\in\mathbb{C}_r^{m×n}$，$A^HA$的特征值满足$\lambda_1\geq\lambda_2\geq···\geq\lambda_r>0$，$\lambda_{r+1}=\lambda_{r+2}=···=\lambda_n=0$，称$\sigma_i=\sqrt{\lambda_i}$为$A$的奇异值。称$\sigma_i,i=1,···,r$为$A$的正奇异值。

11. 由引理1.1.3可知，要计算一个矩阵$A$的正奇异值，只需要计算$A^HA$或者$AA^H$即可。

12. ==思考== ：复方阵的特征值和奇异值有何异同？

    同一矩阵的特征值可能与其奇异值完全相同，可能部分相同，也可能完全不同。

13. ==命题1.1.1== 设$\sigma_1\geq\sigma_2\geq···\geq \sigma_n$和$|\lambda_1|\geq|\lambda_2|\geq···\geq|\lambda_n|$分别是$n$阶正规矩阵$A$（$A^HA=AA^H$）的奇异值和特征值，则$\sigma_i=|\lambda_i|,i=1,···,n$。也就是说正规矩阵的奇异值依顺序等于其特征值的绝对值。

14. ==命题1.1.2== 设$\sigma_i$和$\lambda_i(i=1,···,n)$分别是$n$阶复方阵$A$的奇异值和特征值，则
    $$
    \sum_{i=1}^n|\lambda_i|^2\leq\sum_{i=1}^n|\sigma_i|^2
    $$
    该命题的证明借助了Schur引理。

15. ==定理1.1.1（奇异值分解）== 设$A\in\mathbb{C}_r^{m×n}$，则存在酉矩阵$U\in\mathbb{C}^{m×m}$与酉矩阵$V\in\mathbb{C}^{n×n}$使得
    $$
    A=U\begin{bmatrix}
    \Sigma_r & 0 \\
    0 & 0
    \end{bmatrix}V^H=UDV^H
    $$
    其中$\Sigma_r=diag(\sigma_1,···,\sigma_r)$，$\sigma_i(i=1,···,r)$是矩阵$A$的正奇异值。

    > [!CAUTION]
    >
    > 若$A=UDV^H$是奇异值分解，则酉矩阵$U$的列向量为$AA^H$的特征向量，酉矩阵$V$的列向量为$A^HA$的特征向量。

16. 奇异值分解步骤1：

    > [!IMPORTANT]
    >
    > 1、计算$AA^H$的$m$个奇异值，确定对角阵$\Sigma_r$；
    >
    > 2、计算$AA^H$的$m$个标准正交特征向量$\alpha_1,···,\alpha_m$，将其标准化构成酉矩阵$U$；（这一步需要采用诸如施密特正交化的方法进行）
    >
    > 3、计算$A^H\alpha_1,···,A^H\alpha_r$，将其扩为$\mathbb{C}^n$的一组标准正交基，构成矩阵$V$。
    >
    > $A^H\alpha_1,···,A^H\alpha_r$实际上是$A^HA$对应于$\lambda_1,···,\lambda_r$的正交特征向量，其推导步骤如下：
    >
    > $A^HA(A^H\alpha_1)=A^H(AA^H\alpha_1)=A^H(\lambda_1\alpha_1)=\lambda_1(A^H\alpha_1)$。满足特征向量的定义。

17. 奇异值分解步骤2：

    > [!IMPORTANT]
    >
    > 1、计算$A^HA$的$n$个奇异值，确定对角阵$\Sigma_r$；
    >
    > 2、计算$A^HA$的$n$个标准正交特征向量$\beta_1,···,\beta_n$，将其标准话构成酉矩阵$V$；
    >
    > 3、计算$A\beta_1/\sigma_1,···,A\beta_r/\sigma_r$，将其扩展为$\mathbb{C}^m$的一组标准正交基，构成矩阵$U$。
    >
    > $A\beta_1/\sigma_1,···,A\beta_r/\sigma_r$实际上是$AA^H$对应于$\lambda_1,···,\lambda_r$的正交特征向量。

18. ==扩充标准正交基== ：在求解奇异值分解时，对于酉矩阵$U$或者$V$，有时会出现三维空间只有两个正交基向量。此时需要将其扩充为$\mathbb{R}^3$中的标准正交基，用以满足酉矩阵的定义。

19. 尽管矩阵$U$的列向量为矩阵$AA^H$的特征向量，矩阵$V$的列向量为$A^HA$的特征向量，但是酉矩阵$U$和$V$的选择具有一定相关性。可以先确定酉矩阵$U$，再根据$V_1=A^HU_1\Sigma_r^{-1}$确定$V$的前$r$列；或者先确定酉矩阵$V$，根据$U_1=AV_1\Sigma_r^{-1}$确定$U$的前$r$列。

    > [!CAUTION]
    >
    > 上面所有下标为$r$的部分，都要注意：数量上，$r$为非零奇异值的个数，因此涉及到$r$的计算都不要包含$0$特征值对应的特征向量。

20. ==定理1.1.2（极分解）== 任意复方阵$A$必有如下分解：
    $$
    A=GW
    $$
    式中，$G$为半正定Hermite矩阵，$W$是酉矩阵。当矩阵$A$可逆时，$G$是正定Hermite矩阵，此时该极分解式唯一。

    上述分解式称为左极分解。若$A=VG$，其中$V$是酉矩阵，$G$是半正定Hermite矩阵，则称为右极分解。

21. ==应用（最小二乘问题）== ：奇异值分解在矩阵特征值、广义逆矩阵等矩阵分析和计算方面有着重要应用，而且在图像处理、机器学习等领域有着广泛应用。

### 1.1.2 矩阵方程$AXB=D$

1. ==定义1.2.1（相容矩阵）== 已知矩阵$A\in\mathbb{C}_{r_A}^{m×n}$，$B\in\mathbb{C}_{r_B}^{p×q}$，$D\in\mathbb{C}_{r_D}^{m×q}$，待定矩阵$X\in\mathbb{C}_{r_X}^{n×p}$，定义矩阵方程
   $$
   AXB=D
   $$
   若存在矩阵$X\in\mathbb{C}_{r_X}^{n×p}$使矩阵方程$AXB=D$成立，则称方程$AXB=D$**相容**，否则称为**不相容方程**或者**矛盾方程**。

2. 三种矩阵方程的解法：矩阵分解（相抵分解、满秩分解、奇异值分解）；Penrose定理求解；直积和矩阵拉直求解。

3. <font color='red'>**方法一：利用矩阵分解求解**</font>

   ==定理1.2.1== 考察矩阵方程$AXB=D$，如果存在非奇异矩阵（**指行列式不为零的方阵，具有可逆性**）$P_i$和$Q_i$（$i=A,B$）满足式
   $$
   P_AAQ_A=\begin{bmatrix}
   I_r & 0 \\
   0 & 0
   \end{bmatrix},P_BBQ_B=\begin{bmatrix}
   I_r & 0 \\
   0 & 0
   \end{bmatrix}
   $$
   则方程$AXB=D$相容的充要条件是
   $$
   P_ADQ_B=\begin{bmatrix}
   J_{11} & 0 \\
   0 & 0
   \end{bmatrix}
   $$
   此时方程$AXB=D$的解为$X=Q_A\begin{bmatrix}
   J_r & Y_{12} \\
   Y_{21} & Y_{22}
   \end{bmatrix}P_B$。

   式中，$J_{11}\in\mathbb{C}^{r×r}$；$Y_{12},Y_{21},Y_{22}$是适当阶数的任意矩阵。由于$P_A、D、Q_B$是三个已知的矩阵，所以$J_{11}$可以计算出来。

   > [!NOTE]
   >
   > 证明：将$P_AAQ_A=\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix},P_BBQ_B=\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}$代入方程$AXB=D$可得
   > $$
   > P_A^{-1}\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}Q_A^{-1}XP_B^{-1}\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}Q_B^{-1}=D\\
   > \Rightarrow \\
   > \begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}Q_A^{-1}XP_B^{-1}\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}=P_ADQ_B
   > $$
   > 定义$Y=Q_A^{-1}XP_B^{-1}$，$J=P_ADQ_B$并对其分块
   > $$
   > Y=\begin{bmatrix}
   > 	Y_{11} & Y_{12} \\
   > Y_{21} & Y_{22}
   > \end{bmatrix},
   > J=\begin{bmatrix}
   > J_{11} & J_{12} \\
   > J_{21} & J_{22}
   > \end{bmatrix}
   > $$
   > 式中，$Y_{11}$和$J_{11}\in\mathbb{C}^{r×r}$。由此，式（12）简化为
   > $$
   > \begin{bmatrix}
   > Y_{11} & 0 \\
   > 0 & 0
   > \end{bmatrix}=\begin{bmatrix}
   > J_{11} & J_{12} \\
   > J_{21} & J_{22}
   > \end{bmatrix}=P_ADQ_B
   > $$
   > 所以，方程$AXB=D$相容当且仅当
   > $$
   > P_ADQ_B=\begin{bmatrix}
   > J_{11} & 0 \\
   > 0 & 0
   > \end{bmatrix}
   > $$
   > $P_ADQ_B$是三个已知矩阵的相乘，所以可以计算出$J_{11}$。
   >
   > 此时注意到$Y=Q_A^{-1}XP_B^{-1}$可得矩阵方程的解为
   > $$
   > X=Q_A\begin{bmatrix}
   > J_{11} & Y_{12} \\
   > Y_{21} & Y_{22}
   > \end{bmatrix}P_B
   > $$
   > 式中$J_{11}\in\mathbb{C}^{r×r}$，$Y_{12},Y_{21}$和$Y_{22}$是适当阶数的任意矩阵。

   > [!TIP]
   >
   > 上面的推导可以帮助理解，从解题的角度，我们只需要记住以下几个步骤：
   >
   > Step 1：$P_ADQ_B=\begin{bmatrix}
   > J_{11} & 0 \\
   > 0 & 0
   > \end{bmatrix}$（相容性检查）；
   >
   > Step 2：方程的解$X=Q_A\begin{bmatrix}
   > J_{11} & Y_{12} \\
   > Y_{21} & Y_{22}
   > \end{bmatrix}P_B$。

   > [!CAUTION]
   >
   > 若$r_A=r_B=r$，由定理1.2.1知，方程$AXB=D$相容的必要条件为$r_D\leq r$。也就是说，从$r_D\leq r$可以推导出方程$AXB=D$相容。（因为$D$是由三个矩阵相乘得来，他的秩不可能高于原矩阵的秩。相等的时候分解唯一。）

4. ==推论1.2.1== 考察矩阵方程$AXA=A$，存在非奇异矩阵$P$和$Q$使得$PAQ=\begin{bmatrix}
   I_r & 0 \\
   0 & 0
   \end{bmatrix}$，则$G$是该方程解的充要条件是
   $$
   G=Q\begin{bmatrix}
   I_r & K \\
   L & M
   \end{bmatrix}P
   $$
   式中，$K,L,M$是适当阶数的任意矩阵，当$A$为可逆矩阵时，$G$唯一。

   > [!CAUTION]
   >
   > 推论1.2.1表明矩阵方程$AXA=A$的解一定存在$\Rightarrow$任一矩阵的减号逆一定存在。（可以出判断题）

5. ==例1.2.1== 求解矩阵方程$AXA=A$，其中$A=\begin{bmatrix}
   1&0&1 \\
   0&1&1
   \end{bmatrix}$。`//TODO`

   > [!IMPORTANT]
   >
   > 这一题目的解题过程和解析如下：
   >
   > 解：$\begin{bmatrix}
   > A&I \\
   > I&*
   > \end{bmatrix}=\begin{bmatrix}
   > 1&0&1&1&0 \\
   > 0&1&1&0&1 \\ 1&0&0&*&* \\ 0&1&0&*&* \\ 0&0&1&*&*
   > \end{bmatrix}\Rightarrow(初等变换)\begin{bmatrix}
   > 1&0&0&1&0 \\
   > 0&1&0&0&1 \\ 1&0&-1&*&* \\ 0&1&-1&*&* \\ 0&0&1&*&*
   > \end{bmatrix}$，其中$*$为任意复数。
   >
   > 初等变换包括行变换和列变换，对于可逆的初等变换矩阵，左乘为行变换，右乘为列变换。
   >
   > 由此解得$P=\begin{bmatrix}
   > 1&0 \\
   > 0&1 
   > \end{bmatrix}$，$Q=\begin{bmatrix}
   > 1&0&-1 \\
   > 0&1&-1 \\ 0&0&1
   > \end{bmatrix}$。故方程$AXA=A$的解可以表示为
   > $$
   > G=Q\begin{bmatrix}
   > 1&0 \\
   > 0&1 \\ a&b
   > \end{bmatrix}P=\begin{bmatrix}
   > 1-a&-b \\
   > -a&1-b \\ a&b
   > \end{bmatrix}
   > $$
   > 式中$a,b$为任意数。

   > [!CAUTION]
   >
   > 推论1.2.1中可逆矩阵$P$和$Q$可通过初等变换求解，即
   > $$
   > \begin{bmatrix}
   > P&0 \\
   > 0&I 
   > \end{bmatrix}
   > \begin{bmatrix}
   > A&I \\
   > I&* 
   > \end{bmatrix}
   > \begin{bmatrix}
   > Q&0 \\
   > 0&I 
   > \end{bmatrix}=
   > \begin{bmatrix}
   > PAQ&P \\
   > Q&* 
   > \end{bmatrix}
   > $$
   > **为什么？**（初等变换矩阵均为可逆矩阵，而且根据推论1.2.1，$PAQ=\begin{bmatrix}
   > I_r & 0 \\
   > 0 & 0
   > \end{bmatrix}$，可以通过初等变换求解可逆矩阵$P$和$Q$）

6. ==定理1.2.2== 考察矩阵方程$AXB=D$，若矩阵$A$和$B$有奇异值分解
   $$
   A=U_A\begin{bmatrix}
   \Sigma_r & 0 \\
   0 & 0
   \end{bmatrix}V_A^H,B=U_B\begin{bmatrix}
   \Lambda_r & 0 \\
   0 & 0
   \end{bmatrix}V_B^H
   $$
   则方程$AXB=D$相容的充分必要条件是
   $$
   U_A^HDV_B=\begin{bmatrix}
   J_{11} & 0 \\
   0 & 0
   \end{bmatrix}
   $$
   此时方程的解为
   $$
   X=V_A\begin{bmatrix}
   \Sigma_r^{-1}J_{11}\Lambda_r^{-1} & Y_{12} \\
   Y_{21} & Y_{22}
   \end{bmatrix}U_B^H
   $$

7. ==推论1.2.2== 考察矩阵方程$AXA=A$，若矩阵$A$有奇异值分解$A=U\begin{bmatrix}
   \Sigma_r^{-1} & 0 \\
   0 & 0
   \end{bmatrix}V^H$，则矩阵$G$是该方程解的充要条件是
   $$
   G=V\begin{bmatrix}
   \Sigma_r^{-1} & K \\
   L & M
   \end{bmatrix}U^H
   $$
   式中$\Sigma_r$是$r$阶对角矩阵，其对角线元素是矩阵$A$的$r$个正奇异值；$K,L,M$是适当阶数的任意矩阵。当$A$为可逆矩阵时，$G$唯一。

8. <font color='red'>**方法二：利用Penrose定理求解**</font>

   ==定理1.2.3（Penrose定理）== 矩阵方程$AXB=D$相容的充要条件为$AA^{(1)}DB^{(1)}B=D$，其中$A^{(1)}\in A\{1\}$，$B^{(1)}\in B\{1\}$；此时方程$AXB=D$的通解为
   $$
   X=A^{(1)}DB^{(1)}+Y-A^{(1)}AYBB^{(1)}
   $$
   式中，$Y$为任一$n×p$矩阵。
   $$
   AA^{(1)}A=A,BB^{(1)}B=B
   $$
   找到一个减号逆，遍历矩阵$Y$，即可构造出减号逆的集合$X=A^{(1)}DB^{(1)}+Y-A^{(1)}AYBB^{(1)}$。

   > [!NOTE]
   >
   > Penrose定理证明：
   >
   > 1、充分性（条件成立推结论）：设$AA^{(1)}DB^{(1)}B=D$成立，则参考矩阵方程$AXB=D$的形式，有$X=A^{(1)}DB^{(1)}$使得矩阵方程$AXB=D$相容，即充分性成立；
   >
   > 2、必要性（结论成立推条件）：假设矩阵方程$AXB=D$有解，不妨设$X$是其中一个解。根据推论1.2.1，任何一个矩阵都有其减号逆，则对于矩阵$A$、$B$来说，有$AA^{(1)}A=A$，$BB^{(1)}B=B$，代入原方程可得
   > $$
   > AA^{(1)}AXBB^{(1)}B=D
   > $$
   > 因为$AXB=D$，所以有$AA^{(1)}DB^{(1)}B=D$。必要性成立。
   >
   > 通解的证明：
   >
   > 1、充分性（条件成立推结论）：假设条件成立，即通解为定理中的形式，则将原公式左侧代入，得
   > $$
   > AXB=AA^{(1)}DB^{(1)}B+AYB-AA^{(1)}AYBB^{(1)}B=\\ AA^{(1)}DB^{(1)}B+AYB-AYB=AA^{(1)}DB^{(1)}B=D
   > $$
   > 通解形式的充分性成立。
   >
   > 2、必要性（结论成立推条件）：
   >
   > 假设$X$是矩阵方程$AXB=D$的任意解，则设
   > $$
   > X=A^{(1)}DB^{(1)}+X-A^{(1)}DB^{(1)}
   > $$
   > 将$AXB=D$代入上述公式，有
   > $$
   > X=A^{(1)}DB^{(1)}+X-A^{(1)}AXBB^{(1)}
   > $$
   > 符合通解形式。因此必要性成立。

9. ==推论1.2.3== 非齐次线性方程组$Ax=b$相容的充要条件是$AA^{(1)}b=b$，其通解为
   $$
   x=A^{(1)}b+(I-A^{(1)}A)y, \forall y\in\mathbb{C}^n
   $$

10. ==推论1.2.4== 齐次线性方程组$Ax=0$的通解为
    $$
    x=(I-A^{(1)}A)y, \forall y\in\mathbb{C}^n
    $$

11. ==推论1.2.5== 设$A\in\mathbb{C}^{m×n}$，则
    $$
    A\{1\}=A^{(1)}+Z-A^{(1)}AZAA^{(1)}
    $$
    其中，$Z$为任一$n×m$矩阵。

    关于上述三个推论的例题：==例1.2.2==

12. <font color='red'>**方法三：利用直积和矩阵拉直求解**</font>

    ==定义1.2.2（Kronecker积）== 设$A=(a_{ij})\in\mathbb{C}^{m×n}$，$B=(b_{ij})\in\mathbb{C}^{p×q}$，称如下<u>**分块矩阵**</u>
    $$
    A\otimes B=\begin{bmatrix}
    a_{11}B&a_{12}B&···&a_{1n}B \\
    a_{21}B&a_{22}B&···&a_{2n}B \\ \vdots&\vdots&&\vdots \\ a_{m1}B&a_{m2}B&···&a_{mn}B
    \end{bmatrix}\in\mathbb{C}^{mp×nq}
    $$
    为$A$与$B$的Kronecker积（或直积、张量积），简记为$A\otimes B=(a_{ij}B)$。

    辅助理解的例题：==例1.2.3==

13. ==命题1.2.1（直积的性质）== 矩阵的Kronecker积具有以下性质：

    > [!IMPORTANT]
    >
    > （1）对任意复数 $k$，$(kA) \otimes B = A \otimes (kB) = k(A \otimes B)$;
    >
    > （2）$(A \otimes B) \otimes C = A \otimes (B \otimes C)$;
    >
    > （3）$A \otimes (B + C) = A \otimes B + A \otimes C$;
    >
    > （4）$(A \otimes B)^H = A^H \otimes B^H$;
    >
    > （5）若矩阵 $A$ 和 $C$，矩阵 $B$ 和 $D$ 均可相乘，则 $(A \otimes B)(C \otimes D) = (AC) \otimes (BD)$;
    >
    > （6）$\operatorname{rank}(A \otimes B) = \operatorname{rank}(A)\operatorname{rank}(B)$;
    >
    > （7）$(A \otimes B)^+ = A^+ \otimes B^+$.

14. ==定理1.2.4== 设$f(x,y)=\sum_{i,j=0}^{k}c_{ij}x^iy^j$是变量$x,y$的复二元多项式，对任意矩阵$A\in\mathbb{C}^{m×m}$和$B\in\mathbb{C}^{n×n}$，定义矩阵
    $$
    f(A,B)=\sum_{i,j=0}^{k}c_{ij}(A^i\otimes  B^j)\in\mathbb{C}^{mn×mn}
    $$

    式中，$A^0=I_m$，$B^0=I_n$。若$A$和$B$的特征值分别为$\lambda_1,...,\lambda_m$和$\mu_1,...,\mu_n$，则$f(A,B)$的特征值为$f(\lambda_i,\mu_j)$，$i=1,...,m$，$j=1,...,n$。

15. ==推论1.2.6== 设$A\in\mathbb{C}^{m×m}$，$B\in\mathbb{C}^{n×n}$，其特征值分别为$\lambda_1,...,\lambda_m$和$\mu_1,...,\mu_n$，则有

    （1）$A\otimes B$的特征值为$\lambda_i\mu_j$，$i=1,...,m$，$j=1,...,n$。

    （2）$A\otimes I_n+I_m\otimes B$的特征值为$\lambda_i+\mu_j$，$i=1,...,m$，$j=1,...,n$。

16. ==定义1.2.3（列拉直）== 设$A=(a_{ij})\in\mathbb{C}^{m×n}$，并记$\boldsymbol{a_i}=[a_{1i},a_{2i},...,a_{mi}]^T$（$\boldsymbol{a_i}$此时是一个列向量），$i=1,2,...,n$。令
    $$
    vec(A)=\begin{bmatrix}
    \boldsymbol{a_1} \\
    \boldsymbol{a_2} \\
    \vdots \\
    \boldsymbol{a_n}
    \end{bmatrix}
    $$
    称为矩阵$A$的列拉直（或列展开）。（$vec(A)$是一个列向量，说白了就是将矩阵$A$的每一列前后接起来形成一个新的列向量。）

17. ==定理1.2.5== 设$A\in\mathbb{C}^{m×n}$，$B\in\mathbb{C}^{n×p}$和$C\in\mathbb{C}^{p×q}$，则
    $$
    vec(ABC)=(C^T\otimes A)vec(B)
    $$

18. ==定理1.2.6== 矩阵方程$AXB=D$相容的充要条件为
    $$
    rank(B^T\otimes A,vec(D))=rank(B^T\otimes A)
    $$
    定理1.2.3和定理1.2.6分别给出了矩阵方程$AXB=D$相容的充要条件。显然，这两个条件应该是互为等价的。

    这里可以理解为增广矩阵的秩等于系数矩阵的秩。这里是针对矩阵方程的推广。

19. ==应用== Lyapunov方程，定义为$AX+XA^H=-Q$，该方程在控制理论、通讯和动力系统中起着非常重要的作用.。我们常根据Lyapunov矩阵方程的解来检测系统的稳定性、可控性和可观测性等问题。

20. ==定理1.2.7== 设$A \in \mathbb{C}^{n \times n}$，$Q \in \mathbb{C}^{n \times n}$，若矩阵$A$的所有特征值均具有负实部，则矩阵方程$AX + XA^H = -Q$ 有唯一解，且解$X \in \mathbb{C}^{n \times n}$ 可表示为
    $$
    X=\int_{0}^{+\infty}e^{At}Qe^{A^H t}dt
    $$

    进一步，若 \(Q\) 是（半正定、正定）Hermite矩阵，则解$X$也是（半正定、正定）Hermite矩阵.

### 1.1.3 减号逆

1. 推论1.2.1和推论1.2.2给出了矩阵减号逆的存在性结论，并指出任一矩阵的减号逆必存在。定理1.2.3也提供了矩阵减号逆的通解。此外，我们还可以利用**满秩分解**给出减号逆的一个特解。特别地，若$A\in\mathbb{C}_n^{m×n}$，则$A$的左逆$(A^HA)^{-1}A^H\in A\{1\}$；若$A\in\mathbb{C}_m^{m×n}$，则$A$的右逆$A^H(AA^H)^{-1}\in A\{1\}$。

   （矩阵列满秩有左逆，行满秩有右逆。）

   > [!NOTE]
   >
   > 为什么呢？因为矩阵列满秩的时候，有$m\geq n$，则$A^TA_{n}$为可逆矩阵，有$(A^TA)^{-1}A^TA=I$。所以$A$列满秩时，左逆为$(A^TA)^{-1}A^T$。行满秩右逆同理。

2. ==注1：== 减号逆不唯一。若$A\in\mathbb{C}^{m×n}_r$，则其减号逆阶数为$n×m$（指矩阵的尺寸）。

3. ==例1.3.1== 求矩阵$A=\begin{bmatrix}
   1&0&1 \\
   0&1&1
   \end{bmatrix}$的一个减号逆。

   > [!IMPORTANT]
   >
   > 可以看到矩阵$A$为行满秩矩阵，则矩阵有右逆$A^T(AA^T)^{-1}\in A\{1\}$，为$A$的一个减号逆。
   >
   > 为$A=\begin{bmatrix}
   > \frac{2}{3}&-\frac{1}{3} \\ 
   > -\frac{1}{3}&\frac{2}{3} \\ \frac{1}{3}&\frac{1}{3}
   > \end{bmatrix}$。计算二阶矩阵逆的方法：主对调负反号，除以行列式。

4. ==注2：== 尽管利用奇异值分解或满秩分解求解矩阵的减号逆表达比较简洁，但利用相抵分解求解的计算量相对较小。

5. ==命题1.3.1（减号逆的性质）==

   > [!IMPORTANT]
   >
   > （1）$rank(A)\leq rank(A^-)$；
   >
   > （2）$AA^-$和$A^-A$均为幂等矩阵，且$rank(AA^-)=rank(A^-A)=rank(A)$；
   >
   > ​	$rank(A)\geq rank(AA^-)\geq rank(AA^-A)=rank(A)$；
   >
   > （3）若$B=P^{-1}AQ^{-1}$，则$QA^-P\in B\{1\}$。

6. ==【应用】== 考察非齐次线性方程组$A\boldsymbol{x}=\boldsymbol{b}$。式中，$A\in\mathbb{C}^{m×n}$和$\boldsymbol{b}\in\mathbb{C}^m$给定，$\boldsymbol{x}\in\mathbb{C}^n$为待定向量。

   首先讨论$A^{(1)}$与方程组$A\boldsymbol{x}=\boldsymbol{b}$相容性的关系。

7. ==命题1.3.2（线性方程组相容条件）== 设$A\in\mathbb{C}^{m×n}$，$\boldsymbol{b}\in\mathbb{C}^m$，考察非齐次线性方程组$A\boldsymbol{x}=\boldsymbol{b}$，则以下表达式等价：

   > [!IMPORTANT]
   >
   > （1）$rank(A)=rank(A,\boldsymbol{b})$
   >
   > （2）$\boldsymbol{b}\in R(A)$
   >
   > （3）$AA^{(1)}\boldsymbol{b}=\boldsymbol{b}$（Penrose定理）

8. ==定理1.3.1（减号逆与方程解的关系）== 设$A\in\mathbb{C}^{m×n}$，则$G\in A\{1\}$的充要条件是$\boldsymbol{x}=G\boldsymbol{b}$是相容方程组$A\boldsymbol{x}=\boldsymbol{b}$的解。`//TODO` 找一下证明

9. ==推论1.3.1（减号逆与方程解的关系）== 设$A\in \mathbb{C}^{m×n}$，非零向量$\boldsymbol{b}\in\mathbb{C}^m$，则相容方程组$A\boldsymbol{x}=\boldsymbol{b}$的解一定可以用$\boldsymbol{x}=G\boldsymbol{b}$表示，其中$G\in A\{1\}$。

10. ==例1.3.2== 考察$\boldsymbol{x}=A^{(1)}\boldsymbol{b}$与相容方程组$A\boldsymbol{x}=\boldsymbol{b}$解的关系，其中$A=\begin{bmatrix}
    1&0&0 \\
    0&1&0
    \end{bmatrix}$，$\boldsymbol{b}=\begin{bmatrix}
    1 \\
    0
    \end{bmatrix}$。

    > [!IMPORTANT]
    >
    > Step（1）由于$\boldsymbol{b}\in R(A)$，故方程组$A\boldsymbol{x}=\boldsymbol{b}$相容。（命题1.3.2，线性方程组相容条件）
    >
    > Step（2）通过求解线性方程组（常规求法），可得该方程组的解为$\boldsymbol{x}=\begin{bmatrix}
    > 1 \\
    > 0 \\ c
    > \end{bmatrix}$，其中$\forall c\in \mathbb{C}$。
    >
    > Step（3）由减号逆的定义（直接把$\boldsymbol x$的四个元素设出来计算即可），可以求解得到$A^{(1)}=\begin{bmatrix}
    > 1&0 \\
    > 0&1 \\ a&b
    > \end{bmatrix}$，$\forall a,b\in \mathbb{C}$。则有$A^{(1)}\boldsymbol{b}=\begin{bmatrix}
    > 1 \\
    > 0 \\ a
    > \end{bmatrix}$，$\forall a\in\mathbb{C}$。可以看到，该解的形式和$\boldsymbol{x}$一模一样。那么$\boldsymbol{x}=A^{(1)}\boldsymbol b$。
    >
    > 另外，若$\boldsymbol b=0$，则不论$A^{(1)}$是任何矩阵，方程的解都是$0$。此时，齐次方程组的解为$\boldsymbol x=[0,0,c]^T,\forall c\in\mathbb{C}$。相容方程组的解不能用$\boldsymbol x=A^{(1)}\boldsymbol b$来表示。

11. ==定理1.3.2（相容线性方程组唯一解条件）== 设$A\in\mathbb{C}^{m×n}$，$A^{(1)}A=I_n$当且仅当矩阵$A$是列满秩的。此时，相容方程组$A\boldsymbol{x}=\boldsymbol{b}$有唯一解。（其实就是以前学过的列满秩系数矩阵构成的线性方程组具有唯一解）

    > [!NOTE]
    >
    > 证明：此处证的是当且仅当的部分。
    >
    > 充分性（条件成立推结论）：由命题1.3.1减号逆的性质第二条知，$A^{(1)}A$是幂等矩阵。根据幂等矩阵的性质，它的特征值只能是$0$或$1$。又因为$A$是列满秩矩阵，所以特征值只能全部是$1$。这样，$A^{(1)}A$相似于单位矩阵，即$A^{(1)}A=I_n$。
    >
    > 必要性（结论成立推条件）：设$A^{(1)}A=I_n$，矩阵$A$有左逆$A^{(1)}$，根据**定理3.2.3（教材P99定理+证明）**，矩阵$A$有左逆的充要条件是矩阵$A$是列满秩矩阵。

12. ==注3：== 当矩阵$A$是列满秩矩阵时，它的任一减号逆$A^{(1)}$都是$A$的左逆。当矩阵$A$是行满秩矩阵时，它的任一减号逆$A^{(1)}$都是$A$的右逆。

13. ==例1.3.3== 求向量$\boldsymbol b\in\mathbb{C}^m$在$R(A)$上的正交投影，其中$A\in\mathbb{C}^{m×n}$。（教材P248 例5.3.3）

    这里几个概念要先捋清楚。

    ==正交投影矩阵：==回顾定义3.7.2（幂等矩阵），幂等矩阵也称为投影矩阵。Hermite幂等矩阵称为正交投影矩阵。

    在教材P131的例3.7.2，矩阵$A$是列满秩矩阵，可以设$P=A(A^HA)^{-1}A^H$，$P^2=P$，为幂等矩阵。

    但是在本例中，并没有说明矩阵$A$是不是列满秩矩阵，因此$(A^HA)^{-1}$不一定存在。

    我们把$(A^HA)^{-1}$更换为减号逆$(A^HA)^-$，并令$P=A^H(A^HA)^-A$。

14. ==定理1.3.3== 对任意矩阵$A\in\mathbb{C}^{m×n}$，$A(A^HA)^{-}A^H$是正交投影矩阵。

15. ==注4：== 对一般的矩阵$A$，$(AA^H)^-=(A^H)^-A^-$。

16. ==推论1.3.2== 对任意矩阵$A\in\mathbb{C}^{m×n}$，$A(A^HA)^-A^HA=A$。（证明见教材P249 推论5.3.2）。

    > [!TIP]
    >
    > 另一种证明方法：
    >
    > 设 $G = (A^H A)^-$，已知 $A^H A G A^H A = A^H A$。
    >
    > 考虑 $A - A G A^H A$，左乘 $A^H$：
    >
    > $$
    > A^H (A - A G A^H A) = A^H A - A^H A G A^H A = 0
    > $$
    >
    > 所以 $R(A - A G A^H A) \subseteq N(A^H)$。
    >
    > 但显然又有 $R(A - A G A^H A) \subseteq R(A)$。
    >
    > 而 $R(A) \cap N(A^H) = \{0\}$（这是基本空间的正交补关系），因此：
    >
    > $$
    > A - A G A^H A = 0 \implies A(A^H A)^- A^H A = A
    > $$

17. ==例1.3.4（最小二乘问题，教材P249 例5.4.4（书上序号打错了，应该是5.3.4））== 考察线性方程组$A\boldsymbol x=\boldsymbol b$，其中，矩阵$A\in\mathbb{C}^{m×n}$和向量$\boldsymbol b\in\mathbb{C}^m$给定，向量$\boldsymbol x\in\mathbb{C}^n$待定。求向量$\boldsymbol x$使得$||A\boldsymbol x-\boldsymbol b||_2$最小。

    最小二乘问题是指$b$在$A$的列空间上的正交投影。

### 1.1.4 极小范数广义逆

1. ==定理1.4.1== （教材P250 定理5.4.1，无证明）设矩阵$A\in\mathbb{C}^{m×n}$的奇异值分解为
   $$
   A=U\begin{bmatrix}
   \Sigma_r & 0 \\
   0 & 0
   \end{bmatrix}V^H
   $$
   则$G\in A\{1,4\}$的充要条件是
   $$
   G=V\begin{bmatrix}
   \Sigma_r^{-1} & K \\
   0 & M
   \end{bmatrix}U^H
   $$
   式中，$K$和$M$为适当阶数的任意矩阵。

2. ==注1：== 若$A\in\mathbb{C}^{m×n}_n$，则$(A^HA)^{-1}A^H\in A\{1,4\}$；若$A\in\mathbb{C}^{m×n}_m$，则$A^H(AA^H)^{-1}\in A\{1,4\}$。

3. 【思考】能否利用Penrose定理求解矩阵的极小范数广义逆？第一个条件可以使用Penrose定理，但是第四个条件$(XA)^H=XA$很难使用Penrose定理求解。为了解决该问题，引入定理1.4.2。

4. ==定理1.4.2== （教材P250 定理5.4.2，证明仍然从充分性和必要性两个方向来证明）设$A\in\mathbb{C}^{m×n}$，则
   $$
   A\{1,4\}=\{X\in\mathbb{C}^{n×m}|XA=A_m^-A\}
   $$
   该定理是指，$A\{1,4\}$满足方程$XA=A_m^-A$，$X$是极小范数广义逆的集合中的元素，可以把$A_m^-$看做一个特解，用该特解去检验其他的$X$是否属于集合$A\{1,4\}$。

5. ==定理1.4.3== （教材P250 定理5.4.3）设$A_m^-\in A\{1,4\}$，则
   $$
   A\{1,4\}=\{A_m^-+Z(I-AA_m^-)|Z\in\mathbb{C}^{n×m}\}
   $$

   > [!NOTE]
   >
   > 证明：
   >
   > 由Penrose定理可得，方程$XA=A_m^-A$的通解为
   > $$
   > X=A_m^-AA_m^-+Y-YAA_m^-,\forall Y\in\mathbb{C}^{n×m}
   > $$
   > 设$Z$为$\mathbb{C}^{n×m}$中任意矩阵，将$Y=A_m^-+Z$代入上式可得
   > $$
   > X=A_m^-+Z(I-AA_m^-)
   > $$

6. 【应用】当相容矩阵$A\boldsymbol x=\boldsymbol b$有多个解时，往往需要从中挑选一个合适的解。此时我们一般挑选所有解中范数最小的解，或者是其中之一，即如下优化问题
   $$
   \mathop{\mathop{min}\limits_{A\boldsymbol x=\boldsymbol b}}\limits_{AA^{(1)}\boldsymbol b=\boldsymbol b}||\boldsymbol x||_2
   $$
   将满足该优化问题的解称为相容线性方程组$A\boldsymbol x=\boldsymbol b$的极小范数解。

7. ==定理1.4.4== （教材P251 定理5.4.4，注意看一下应用后续的问题建模）矩阵$G\in A\{1,4\}$的充要条件为$\boldsymbol x=G\boldsymbol b$是相容方程组$A\boldsymbol x=\boldsymbol b$的极小范数解。

   该定理和定理5.3.1类似，定理5.3.1指出，设$A\in\mathbb{C}^{m×n}$，则$G\in A\{1\}$的充要条件是$\boldsymbol{x}=G\boldsymbol{b}$是相容方程组$A\boldsymbol{x}=\boldsymbol{b}$的解。

   > [!NOTE]
   >
   > 证明：
   >
   > 必要性（由结论推条件）：在该方向上，设$G\in A\{1,4\}$成立。则首先设$G\in A\{1\}$，那么相容方程$A\boldsymbol x=\boldsymbol b$的通解为
   > $$
   > \boldsymbol x=G\boldsymbol b+(I-GA)\boldsymbol y,\forall y\in\mathbb{C}^n
   > $$
   > 这里注意，因为我们已知矩阵方程是相容的，所以$G$作为$A$的减号逆之一，满足矩阵方程相容的充要条件$AA^{(1)}\boldsymbol b=\boldsymbol b$。即有式$(48)$的通解。
   >
   > 因为矩阵方程是相容的，$\therefore \exist \boldsymbol z \in \mathbb{C}^n$使得$A\boldsymbol z=\boldsymbol b$。
   >
   > 又因为$G$满足第四条，即$(GA)^H=GA$，那么我们求$\boldsymbol x$的二范数，使得$||\boldsymbol x||_2^2$最小。表达式为
   > $$
   > ||\boldsymbol x||_2^2=||G\boldsymbol b+(I-GA)\boldsymbol y||^2_2\\
   > =||G\boldsymbol b||_2^2+||G\boldsymbol b+(I-GA)\boldsymbol y||^2_2+2(G\boldsymbol b)·[(I-GA)\boldsymbol y]
   > $$
   > 因为$G\boldsymbol b$为一个列向量，所以它和后面一项做点乘可以写作
   > $$
   > (G\boldsymbol b)^H(I-GA)\boldsymbol y=(GA\boldsymbol z)^H(I-GA)\boldsymbol y\\
   > =\boldsymbol z^HGA(I-GA)\boldsymbol y\\
   > =\boldsymbol z^H(GA-GAGA)\boldsymbol y=0
   > $$
   > 因为$G\in A\{1,4\}$，所以$AGA=A$。此时有
   > $$
   > ||\boldsymbol x||_2^2=||G\boldsymbol b||_2^2+||G\boldsymbol b+(I-GA)\boldsymbol y||^2_2\geq ||G\boldsymbol b||_2^2
   > $$
   > 充分性（由条件推结论）：在该方向上，假设已知$\boldsymbol x=G\boldsymbol b$为相容方程组的极小范数解。根据==定理5.3.1==，$G\in A\{1\}$。
   >
   > 根据范数不等式，有
   > $$
   > ||\boldsymbol x||_2^2=||G\boldsymbol b+(I-GA)\boldsymbol y||^2_2\\
   > \leq ||G\boldsymbol b||_2^2+||(I-GA)\boldsymbol y||^2_2\\
   > \leq ||G\boldsymbol b||_2^2
   > $$
   >
   > $\because$等式严格成立，$\therefore$需要保证$(G\boldsymbol b)·[(I-GA)\boldsymbol y]=0,(I-GA)\boldsymbol y=0$。那么有
   > $$
   > (G\boldsymbol b)^H[(I-GA)\boldsymbol y]=(GA\boldsymbol z)^H[(I-GA)\boldsymbol y]\\
   > =\boldsymbol z^H(GA)^H(I-GA)\boldsymbol y\\
   > =z^H[(GA)^H-(GA)^HGA]\boldsymbol y=0
   > $$
   > 当$(GA)^H=GA$时，上式成立。所以$G\in A\{4\}$。结合第一条，我们有$G\in A\{1,4\}$。

   该定理在解题中的应用：如果我们通过其他的方法获得了$G\in A\{1,4\}$，通过该定理可以直接得出$\boldsymbol x=G\boldsymbol b$为相容方程组$A\boldsymbol x=\boldsymbol b$的极小范数解。

8. ==推论1.4.1== （教材P252 推论5.4.1）相容方程$A\boldsymbol x=\boldsymbol b$的极小范数解是唯一的。

   > [!NOTE]
   >
   > 证明：设$G_1\boldsymbol b$和$G_2\boldsymbol b$都是方程组$A\boldsymbol x=\boldsymbol b$的极小范数解，则$G_1,G_2\in A\{1,4\}$。又得知线性方程组$A\boldsymbol x=\boldsymbol b$是相容的，必存在向量$\boldsymbol z\in\mathbb{C}^n$使得$A\boldsymbol z=\boldsymbol b$。此时有
   > $$
   > G_1\boldsymbol b=G_1A\boldsymbol z,G_2\boldsymbol b=G_2A\boldsymbol z
   > $$
   > 根据定理5.4.2知，$G_1A=G_2A$。因此，$G_1\boldsymbol b=G_2\boldsymbol b$。

9. ==例1.4.1== 求相容方程组$A\boldsymbol x=\boldsymbol b$的极小范数解，其中$A=\begin{bmatrix}
   1&2&0 \\
   2&0&2
   \end{bmatrix}$，$\boldsymbol b=\begin{bmatrix}
   1 \\
   1
   \end{bmatrix}$。

   > [!IMPORTANT]
   >
   > $\because A$行满秩，$\therefore G=A^H(AA^H)^{-1}\in A\{1,4\}$，根据注1。
   >
   > 再根据定理5.4.4，则有$x=G\boldsymbol b$为相容方程组的极小范数解。最终解得$\boldsymbol x=\begin{bmatrix}\frac{1}{3} \\
   > \frac{1}{3}\\ \frac{1}{6}\end{bmatrix}$。
   >
   > 逻辑链：找一个特解$G\in A\{1,4\}$，再通过$\boldsymbol x=G\boldsymbol b$求解。
   >
   > 另一种解法，奇异值分解（见教材P252、P253解法二，根据定理1.4.1求解）。

10. ==例1.4.2== 求向量$\boldsymbol b\in\mathbb{C}^m$在$R(A^H)$上的正交投影，其中$A\in\mathbb{C}^{m×n}$。

    > [!IMPORTANT]
    >
    > 根据定理1.3.3，我们取正交投影矩阵$P=A^H(AA^H)^-A$，那么$P\boldsymbol b$即为向量$\boldsymbol b$在$R(A^H)$上的正交投影。$P\boldsymbol b=A_m^-A\boldsymbol b$。
    >
    > 这道题要看教材P253 ==例5.4.2==，通过奇异值分解来求解和证明的。

### 1.1.5 最小二乘广义逆

最小二乘广义逆是指满足Penrose条件第1条和第3条的逆。

1. ==定理1.5.1== 设$A\in\mathbb{C}^{m×n}_r$的奇异值分解为
   $$
   A=U\begin{bmatrix}
   \Sigma_r & 0 \\
   0 & 0
   \end{bmatrix}V^H
   $$
   则$G\in A\{1,3\}$的充要条件是
   $$
   G=V\begin{bmatrix}
   \Sigma_r^{-1} & 0 \\
   L & M
   \end{bmatrix}U^H
   $$
   式中，$L$和$M$为适当阶数的任意矩阵。

2. ==注1：== 若$A\in\mathbb{C}^{m×n}_n$，则$(A^HA)^{-1}A^H\in A\{1,3\}$；若$A\in\mathbb{C}^{m×n}_m$，则$A^H(AA^H)^{-1}\in A\{1,3\}$。

3. ==定理1.5.2== （教材P254 定理5.4.2）设$A\in\mathbb{C}^{m×n}$，则
   $$
   A\{1,3\}=\{X\in\mathbb{C}^{n×m}|AX=AA_l^-\}
   $$

4. ==定理1.5.3== （教材P254 定理5.4.3）设$A_l^-\in A\{1,3\}$，则
   $$
   A\{1,3\}=\{A_l^-+(I-A_l^-A)Z|Z\in\mathbb{C}^{n×m}\}
   $$
   证明通过Penrose定理进行。该定理是说，当我们已知一个最小二乘广义逆时，我们可以使用该最小二乘广义逆去写出最小二乘广义逆的集合。

5. 【应用，线性最小二乘问题】矩阵的最小二乘广义逆与最小二乘问题有密切关系。考察线性方程组$A\boldsymbol x=\boldsymbol b$，式中，$A\in\mathbb{C}^{m×n}$和$\boldsymbol b\in\mathbb{C}^m$给定，$\boldsymbol x\in\mathbb{C}^n$待定。若方程组不相容，则他没有通常意义上的解。退而求其次我们求使得残量范数（向量2范数）最小的解，即求向量$\boldsymbol x$使得$||A\boldsymbol x-\boldsymbol b||_2$最小。通常将这一优化问题称为线性最小二乘问题，向量$\boldsymbol x$为不相容方程组的最小二乘解。

6. ==定理1.5.4== （证明见教材P255 定理5.4.4）对任意矩阵$A\in\mathbb{C}^{m×n}$，$AA_l^-$是正交投影矩阵。

7. ==推论1.5.1== 向量$\boldsymbol b$在$R(A)$的正交投影为$AA_l^-\boldsymbol b$，其中$A\in\mathbb{C}^{m×n}$和$\boldsymbol b\in\mathbb{C}^m$。

8. ==推论1.5.2== 给定$A\in\mathbb{C}^{m×n}$和$\boldsymbol b\in\mathbb{C}^m$，有$A(A^HA)^-A^H\boldsymbol b=AA_l^-\boldsymbol b$。

9. ==定理1.5.6== 矩阵$G\in A\{1,3\}$的充要条件为$\boldsymbol x=G\boldsymbol b$是不相容方程$A\boldsymbol x=\boldsymbol b$的最小二乘解。

   > [!NOTE]
   >
   > 证明：
   >
   > 必要性（由结论推条件）：设$G\in A\{1,3\}$，则有$AGA=A$和$(GA)^H=GA$。令$\boldsymbol x=G\boldsymbol b$，则$A\boldsymbol x-\boldsymbol b=AG\boldsymbol b-\boldsymbol b$。
   > $$
   > ||A\boldsymbol x-\boldsymbol b||^2_2=||A\boldsymbol x||^2_2+||\boldsymbol b||^2_2+2(A\boldsymbol x)·\boldsymbol b\\
   > =||AG\boldsymbol b||^2_2+||\boldsymbol b||^2_2+2(AG\boldsymbol b)^H\boldsymbol b\\
   > =
   > $$
   > 这表明$\boldsymbol x=G\boldsymbol b$是不相容方程组$A\boldsymbol x=\boldsymbol b$的最小二乘解。
   >
   > 充分性（由条件推结论）：

10. ==推论1.5.3== 向量$\boldsymbol x$是不相容方程组$A\boldsymbol x=\boldsymbol b$的最小二乘解当且仅当$\boldsymbol x$是相容方程组$A\boldsymbol x=AA_l^-\boldsymbol b$的解，且$A\boldsymbol x=\boldsymbol b$的最小二乘解通式为
    $$
    \boldsymbol x=A_l^-\boldsymbol b+(I-A_l^-A)\boldsymbol y,\forall\boldsymbol y\in\mathbb{C}^n
    $$
    在求解题目时，可以通过其他方法解出一个$A_l^-$，然后通过通式写出$\boldsymbol x$。

11. ==注2：== 方程组$A\boldsymbol x=AA_l^-\boldsymbol b$、$A\boldsymbol x=A(A^HA)^-A^H\boldsymbol b$和$A^HA\boldsymbol x=A^H\boldsymbol b$均是相容同解方程组。

12. ==例1.5.1== 求方程组$A\boldsymbol x=\boldsymbol b$的最小二乘通解式，其中$A=\begin{bmatrix}
    1&2&0 \\
    2&4&0
    \end{bmatrix}$，$\boldsymbol b=\begin{bmatrix}
    1 \\
    1
    \end{bmatrix}$。

    > [!IMPORTANT]
    >
    > 思路：
    >
    > 步骤一：通过定理1.5.1，使用奇异值分解求出一个$G=A_l^-\in A\{1,3\}$；
    >
    > 步骤二：将$G$代入推论1.5.3中的通解式即可。

13. ==定理1.5.5== 不相容方程组$A\boldsymbol x=\boldsymbol b$具有唯一的最小二乘解当且仅当$A$是列满秩矩阵。

    > [!NOTE]
    >
    > 证明：由推论1.5.1知，不相容方程组$A\boldsymbol x=\boldsymbol b$的最小二乘解与相容方程$A\boldsymbol x=AA_l^-\boldsymbol b$等价。再由定理1.3.2可知，相容方程组有唯一解当且仅当$A$列满秩。

14. ==例1.5.2== 求方程组$A\boldsymbol x=\boldsymbol b$的最小二乘解通式，其中$A=\begin{bmatrix}
    1&2 \\
    2&1 \\ 1&1
    \end{bmatrix}$，$\boldsymbol b=\begin{bmatrix}
    1 \\
    0 \\ 0
    \end{bmatrix}$。

    > [!IMPORTANT]
    >
    > 定理1.5.5+注1+定理1.5.6：因为$A$列满秩，所以根据注1，$(A^HA)^-A^H$是最小二乘逆；又因为定理1.5.5，$\boldsymbol x=G\boldsymbol b$是不相容方程组的最小二乘解；又因为定理1.5.6，最小二乘解唯一。所以求出来的最小二乘解唯一。

### 1.1.6 加号逆

加号逆是指满足全部四条Penrose方程的逆。

1. ==定理1.6.1（加号逆存在性定理）== 设\(A\in\mathbb{C}^{m×n}_r\)的奇异值分解为
   $$
   A=U\begin{bmatrix}
   \Sigma_r & 0 \\
   0 & 0
   \end{bmatrix}V^H
   $$
   则
   $$
   G=V\begin{bmatrix}
   \Sigma_r^{-1} & 0 \\
   0 & 0
   \end{bmatrix}U^H
   $$
   是\(A\)的一个加号逆。

   > [!NOTE]
   >
   > **证明：** 将矩阵\(A\)和\(G\)依次代入四个Penrose方程，得
   >
   > $$
   > AGA = U \begin{bmatrix} \Sigma_r & 0 \\ 0 & 0 \end{bmatrix} V^H = A
   > $$
   >
   > $$
   > GAG = V \begin{bmatrix} \Sigma_r^{-1} & 0 \\ 0 & 0 \end{bmatrix} U^H = G
   > $$
   >
   > $$
   > (AG)^H = U \begin{bmatrix} I_r & 0 \\ 0 & 0 \end{bmatrix} U^H = AG
   > $$
   >
   > $$
   > (GA)^H = V \begin{bmatrix} I_r & 0 \\ 0 & 0 \end{bmatrix} V^H = GA
   > $$
   >
   > 因此，矩阵\(G\)是矩阵\(A\)的一个加号逆。

   定理1.6.1利用奇异值分解给出了矩阵$A$的最小二乘广义逆的一个特解，该定理表明矩阵的加号逆一定存在，也表明矩阵的任一广义逆矩阵都存在。

2. ==定理1.6.2（加号逆唯一性定理）== 任给矩阵\(A \in \mathbb{C}^{m \times n}\)，其加号逆\(A^+\)唯一。

   > [!NOTE]
   >
   > 证明：假设\(X\)和\(Y\)均满足四个Penrose方程
   > $$
   > \begin{aligned}
   > X &= XAX \quad (\text{第}2) \\
   >   &= XAYAX = (XA)^H(YA)^HX \quad (\text{第}1、4) \\
   >   &= (AXA)^H Y^H X = A^H Y^H X = YAX \quad (\text{第}1、4) \\
   >   &= YAYAX = Y(AY)^H(AX)^H = YY^H A^H \quad (\text{第}2、3、1) \\
   >   &= YAY \quad (\text{第}3) \\
   >   &= Y \quad (\text{第}2)
   > \end{aligned}
   > $$

3. ==定理1.6.3== （教材P259 定理5.6.3）设矩阵\(A = (a_{ij}) \in \mathbb{C}^{m \times n}_r\)，则
   $$
   A^+ = V_1(\Sigma_r^2)^{-1}V_1^H A^H
   $$
   其中，\(\Sigma_r = \mathrm{diag}(\sigma_1, \cdots, \sigma_r)\)，\(\sigma_1, \cdots, \sigma_r\)是矩阵\(A\)的\(r\)个正奇异值（即非零特征值），\(V_1 \in \mathbb{C}^{n \times r}\)是由\(A^H A\)的属于特征值\(\sigma_1^2, \cdots, \sigma_r^2\)的\(r\)个单位正交特征向量构成。

   该定理直接给出了如何通过公式计算加号逆。

4. ==例1.6.1== 已知矩阵 \(A = \begin{bmatrix} -1 & 0 & 1 \\ 2 & 0 & -2 \end{bmatrix}\)，求 \(A^+\)。

   > [!IMPORTANT]
   >
   > 按照定理1.6.3的方法来直接计算即可。

5. ==推论1.6.1（秩1公式）== 设$A=(a_{ij})\in\mathbb{C}_1^{m×n}$，$A^+=\frac{1}{\sum_{i,j}|a_{i,j}|^2}A^H$。即$A^+=\frac{1}{\sigma_1^2}A^H$。

   证明见教材P259、P260。该推论是由定理1.6.3推理得来。

6. ==定理1.6.4== 设矩阵 \(A \in \mathbb{C}_r^{m \times n}\) 有满秩分解 \(A=BC\)，其中，\(B \in \mathbb{C}_r^{m \times r}\)，\(C \in \mathbb{C}_r^{r \times n}\)，则
   $$
   \begin{aligned}
   A^+ &= C^H (CC^H)^{-1} (B^H B)^{-1} B^H \\
       &= C^H (B^H A C^H)^{-1} B^H
   \end{aligned}
   $$

   > [!NOTE]
   >
   > 证明：将$A=BC$和$G=C^H(CC^H)^{-1}(B^HB)^{-1}B^H$依次代入四个Penrose方程，即可证明。

7. ==推论1.6.2== 若$A\in\mathbb{C}^{m×n}_n$，$A^+=(A^HA)^{-1}A^H$；若$A\in\mathbb{C}^{m×n}_m$，$A^+=A^H(AA^H)^{-1}$。

8. ==例1.6.2== 求 \(A = \begin{bmatrix} 1 & 2 & 1 & 2 \\ 2 & 4 & 2 & 4 \\ 1 & 0 & 0 & 2 \end{bmatrix}\) 的伪逆。

   > [!IMPORTANT]
   >
   > 方法一：奇异值分解
   >
   > 方法二：满秩分解
   >
   > 如何求矩阵的满秩分解呢？
   >
   > 1、将矩阵$A$通过初等行变换化作行阶梯矩阵；
   >
   > 2、取$A$中线性无关的列组成列满秩矩阵$B$；
   >
   > 3、将$A$写作列向量的形式，然后将$A$中所有的列向量在以$B$的列作为基的情况下的坐标写出来；
   >
   > 4、这些坐标写成列向量，组合后即为行满秩矩阵$C$；
   >
   > 5、$A=BC$，分解完毕。

9. ==命题1.6.1（加号逆的性质）== 设 \(A \in \mathbb{C}^{m \times n}\)，则

   > [!IMPORTANT]
   >
   > (1) \((A^+)^+ = A\)；
   >
   > (2) \((A^H)^+ = (A^+)^H\)；
   >
   > (3) \((A^T)^+ = (A^+)^T\)；
   >
   > (4) \(\text{rank}(A^+) = \text{rank}(A)\)；
   >
   > (5) \((A^H A)^+ = A^+(A^H)^+\); \((AA^H)^+ = (A^H)^+A^+\)；
   >
   > (6) \(A^+ = (A^H A)^+ A^H = A^H(AA^H)^+\)。

10. ==注1：== 对于一般的矩阵 \(A\) 和 \(B\)，\((AB)^+ \neq B^+A^+\)；\(A^+A \neq AA^+ \neq I\)。

11. ==推论1.6.3== （教材P262 推论5.6.3，无证明）线性方程组$A\boldsymbol x=\boldsymbol b$相容的充分必要条件为$AA^+\boldsymbol b=\boldsymbol b$，此时通解为
    $$
    \boldsymbol x=A^+\boldsymbol b+(I-A^+A)\boldsymbol y,\forall\boldsymbol y\in\mathbb{C}^n
    $$
    若方程不相容，则上式为其最小二乘通解式。

12. ==定理1.6.5== 向量$\boldsymbol x=A^+\boldsymbol b$既是相容方程$A\boldsymbol x=\boldsymbol b$的唯一极小范数解，也是不相容方程$A\boldsymbol x=\boldsymbol b$的唯一最佳逼近解。

    > [!NOTE]
    >
    > 证明：$\because A^+\in A\{1,4\}$，$\therefore \boldsymbol x=A^+\boldsymbol b$是相容方程的唯一极小范数解（定理1.4.4+推论1.4.1，相容方程的极小范数解是唯一的，且当$G\in A\{1,4\}$的充要条件是$\boldsymbol x=G\boldsymbol b$是相容方程的极小范数解）。
    >
    > 由推论1.5.3可知（向量$\boldsymbol x$是不相容方程组$A\boldsymbol x=\boldsymbol b$的最小二乘解当且仅当$\boldsymbol x$是相容方程组$A\boldsymbol x=AA^+\boldsymbol b$的解，因为$A^+\in A\{1,3\}$，所以此处将原推论中的$A_l^-$更换为$A^+$）。又因为相容方程的唯一极小范数解为$\boldsymbol x=A^+AA^+\boldsymbol b=A^+\boldsymbol b$，所以是唯一的最佳逼近解。

13. ==例1.6.3== 求线性方程 \(Ax=b\) 的极小范数解或最佳逼近解，其中
    $$
    A = \begin{bmatrix}
    1 & 0 & -1 & 1 \\
    0 & 2 & 2 & 2 \\
    -1 & 4 & 5 & 3
    \end{bmatrix}, \quad
    b = \begin{bmatrix}
    4 \\
    -2 \\
    -2
    \end{bmatrix}
    $$

    > [!IMPORTANT]
    >
    > 解题思路：按照奇异值分解（定理1.6.3）或者满秩分解（定理1.6.4）求得$A^+$后，再通过定理1.6.5，计算$\boldsymbol x=A^+\boldsymbol b$即可。

14. ==例1.6.4（最小二乘问题）== 考查线性方程组 \(Ax=b\)，其中，矩阵 \(A \in \mathbb{C}^{m \times n}\) 和向量 \(b \in \mathbb{C}^m\) 给定，向量 \(x \in \mathbb{C}^n\) 待定。求向量 \(x\) 使得 \(\|Ax-b\|_2\) 最小。

    > [!IMPORTANT]
    >
    > \(A^H Ax = A^H b\) 范数最小的解（唯一）
    > $$
    > \begin{gather*}
    > x = (A^H A)^{(1,4)} A^H b \\
    > A = U \begin{bmatrix} \Sigma_r & 0 \\ 0 & 0 \end{bmatrix} V^H \\
    > x = A^{(+)} b \\
    > Ax = AA^{(+)} b
    > \end{gather*}
    > $$

15. ==例1.6.5== （教材P262 例5.6.4）求向量 \(b \in \mathbb{C}^n\) 在 \(R(A)\) 和 \(N(A)\) 上的正交投影，其中 \(A \in \mathbb{C}^{n \times n}\)。

    > [!IMPORTANT]
    >
    > 由于$AA^+$为正交投影矩阵（根据正交投影矩阵和$A^+$的定义得来），且$AA^+\boldsymbol b\in R(A)$，故$Proj_{R(A)}=AA^+\boldsymbol b$。又因为$N(A)\oplus R(A)$是正交直和分解，则$\boldsymbol b\in\mathbb{C}^n$在$N(A)$上的正交投影为
    > $$
    > Proj_{N(A)}\boldsymbol b=[I-A^H(A^H)^+]\boldsymbol b=(I-A^+A)\boldsymbol b
    > $$
    > 上面的等号是通过加号逆的性质第2条实现的。
    > $$
    > \begin{gather*}
    > P = A(A^H A)^- A^H = AA^+ \\
    > Pb = AA^+ b
    > \end{gather*}
    > $$

## 1.2 矩阵微积分及其应用

### 1.2.1 函数矩阵的导数与积分

1. ==定义2.1.1（函数矩阵）== 以变量 $t$ 的函数为元素的矩阵 $A(t)=(a_{ij}(t))_{m\times n}$ 称为函数矩阵。若每个元素 $a_{ij}(t)$ 在 $[a,b]$ 上连续、可微或可积，则称 $A(t)$ 在 $[a,b]$ 上连续、可微或可积。此时定义
   $$
   A'(t)=\frac{d}{dt}A(t)=(a'_{ij}(t))_{m\times n},\quad \int_a^b A(t)dt=\left(\int_a^b a_{ij}(t)dt\right)_{m\times n}
   $$

2. ==命题2.1.1（求导法则）== 设 $A(t),B(t)$ 是适当阶的可微矩阵，$\lambda(t)$ 为可微函数，则
   > [!IMPORTANT]
   > (1) $\frac{d}{dt}(A+B)=\frac{dA}{dt}+\frac{dB}{dt}$；
   >
   > (2) $\frac{d}{dt}(\lambda A)=\frac{d\lambda}{dt}A+\lambda\frac{dA}{dt}$；
   >
   > (3) $\frac{d}{dt}(AB)=\frac{dA}{dt}B+A\frac{dB}{dt}$；
   >
   > (4) 若 $u=f(t)$ 可微，则 $\frac{d}{dt}A(u)=f'(t)\frac{dA}{du}$；
   >
   > (5) 若 $A(t)$ 可逆，则 $\frac{d}{dt}(A^{-1})=-A^{-1}\frac{dA}{dt}A^{-1}$。
   
3. ==命题2.1.2（矩阵指数与三角函数的导数）== 设 $A\in\mathbb{C}^{n\times n}$，则有
   $$
   \frac{d}{dt}e^{At}=Ae^{At}=e^{At}A\\ 
   \quad \frac{d}{dt}\sin At=A\cos At=(\cos At)A\\
   \quad \frac{d}{dt}\cos At=-A\sin At=-(\sin At)A.
   $$
   > [!NOTE]
   > 矩阵指数 $e^{At}$ 的定义可参照常系数线性微分方程组理论，其导数形式与标量情形类似，但注意乘法不可交换，上述公式中 $A$ 与 $e^{At}$ 可交换。
   
4. ==推论2.1.1（积分公式）== 设 $A\in\mathbb{C}^{n\times n}$，则
   $$
   \int_{t_0}^t Ae^{As}ds=e^{At}-e^{At_0}\\
   \quad
   \int_{t_0}^t A\sin As\,ds=\cos At_0-\cos At\\
   \quad
   \int_{t_0}^t A\cos As\,ds=\sin At-\sin At_0.
   $$
   
5. ==命题2.1.3（积分性质）== 设 $A(t),B(t)$ 在 $[a,b]$ 上可积，则
   > [!IMPORTANT]
   > (1) $\int_a^b (A+B)dt=\int_a^b A\,dt+\int_a^b B\,dt$；
   >
   > (2) $\int_a^b \lambda A\,dt=\lambda\int_a^b A\,dt$；
   >
   > (3) 若 $A(t)$ 连续，则 $\frac{d}{dt}\int_a^t A(\tau)d\tau=A(t)$；
   >
   > (4) 若 $A'(t)$ 连续，则 $\int_a^b A'(\tau)d\tau=A(b)-A(a)$。

### 1.2.2 矩阵对矩阵的导数

1. ==定义2.1.2（矩阵对矩阵的导数）== 设 $F(X)=(f_{ij}(X))\in\mathbb{C}^{m\times n}$ 是 $X\in\mathbb{C}^{p\times q}$ 的函数矩阵，每个 $f_{ij}(X)$ 对 $X$ 可微，令
   $$
   \frac{dF(X)}{dX}=\left(\frac{\partial F}{\partial x_{ij}}\right)_{mp\times nq}
   =\begin{bmatrix}
   \frac{\partial F}{\partial x_{11}} & \frac{\partial F}{\partial x_{12}} & \cdots & \frac{\partial F}{\partial x_{1q}} \\
   \frac{\partial F}{\partial x_{21}} & \frac{\partial F}{\partial x_{22}} & \cdots & \frac{\partial F}{\partial x_{2q}} \\
   \vdots & \vdots & \ddots & \vdots \\
   \frac{\partial F}{\partial x_{p1}} & \frac{\partial F}{\partial x_{p2}} & \cdots & \frac{\partial F}{\partial x_{pq}}
   \end{bmatrix},
   \quad\text{其中 }\frac{\partial F}{\partial x_{ij}}=\left(\frac{\partial f_{kl}}{\partial x_{ij}}\right)
   $$
   称 $\frac{dF}{dX}$ 为 $F(X)$ 对 $X$ 的导数。此定义涵盖了向量对向量、向量对矩阵、矩阵对向量和矩阵对矩阵的导数。

2. ==例2.1.2（梯度向量）== 设 $f(x)=f(x_1,\dots,x_n)$ 是 $n$ 元实可微函数，记 $x=[x_1,\dots,x_n]$，则
   $$
   \frac{df}{dx}=\left[\frac{\partial f}{\partial x_1},\frac{\partial f}{\partial x_2},\cdots,\frac{\partial f}{\partial x_n}\right],
   $$
   常称 $\frac{df}{dx}$ 为 $f$ 的梯度向量，记为 $\nabla f$。

3. ==例2.1.3（Jacobian矩阵）== 设 $f(x)=[f_1,\dots,f_m]^T$ 是 $n$ 元实可微向量函数，$x=[x_1,\dots,x_n]$，则
   $$
   \frac{df}{dx}=\begin{bmatrix}
   \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n} \\
   \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n} \\
   \vdots & \vdots & \ddots & \vdots \\
   \frac{\partial f_m}{\partial x_1} & \frac{\partial f_m}{\partial x_2} & \cdots & \frac{\partial f_m}{\partial x_n}
   \end{bmatrix},
   $$
   常称 $\frac{df}{dx}$ 为向量函数 $f$ 的 Jacobian 矩阵。

4. ==补充：Hessian矩阵== 考察 $n$ 元函数 $f(x)$ 在 $x_0$ 的泰勒展开，其二阶项矩阵为
   $$
   H_f(x_0)=\begin{bmatrix}
   \frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1\partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_1\partial x_n} \\
   \frac{\partial^2 f}{\partial x_2\partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \cdots & \frac{\partial^2 f}{\partial x_2\partial x_n} \\
   \vdots & \vdots & \ddots & \vdots \\
   \frac{\partial^2 f}{\partial x_n\partial x_1} & \frac{\partial^2 f}{\partial x_n\partial x_2} & \cdots & \frac{\partial^2 f}{\partial x_n^2}
   \end{bmatrix},
   $$
   该矩阵称为 Hessian 矩阵。此时泰勒展开可写为
   $$
   f(x)=f(x_0)+\nabla f(x_0)(x-x_0)+\frac{1}{2}(x-x_0)^T H_f(x_0)(x-x_0)+\cdots
   $$

### 1.2.3 矩阵微积分的应用

1. ==应用一：最小二乘问题== 考察 $\min_{x\in\mathbb{R}^n}\|Ax-b\|_2^2$，其中 $A\in\mathbb{R}^{m\times n},b\in\mathbb{R}^m$ 给定，$x\in\mathbb{R}^n$ 待定。
   > [!IMPORTANT]
   > 由引理2.1.1（无约束极值的必要条件），极小点需满足 $\frac{d}{dx}\|Ax-b\|_2^2=0$。展开得
   > $$
   > \|Ax-b\|_2^2=x^TA^TAx-2b^TAx+b^Tb,
   > $$
   > 求导：$\frac{d}{dx}(x^TA^TAx)=2A^TAx$，$\frac{d}{dx}(b^TAx)=A^Tb$，故极值条件为
   > $$
   > A^TAx=A^Tb,
   > $$
   > 即最小二乘问题导出的正规方程。这与利用广义逆矩阵理论得到的结论一致。

2. ==应用二：一阶线性常系数微分方程组== 考察
   $$
   \begin{cases}
   \frac{dx(t)}{dt}=Ax(t)+f(t), \\
   x(t_0)=c,
   \end{cases}
   $$
   其中 $A\in\mathbb{C}^{n\times n}$，$f(t)$ 为已知向量函数，$c$ 为初始向量。
   > [!IMPORTANT]
   > 通解公式为
   > $$
   > x(t)=e^{A(t-t_0)}c+\int_{t_0}^t e^{A(t-\tau)}f(\tau)d\tau.
   > $$
   > 若 $f(t)\equiv0$，则解为 $x(t)=e^{A(t-t_0)}c$。

3. ==例2.1.6== 求 $\frac{dx(t)}{dt}=Ax(t)+f(t)$ 的解，其中
   $$
   A=\begin{bmatrix}2&0&0\\1&1&1\\1&-1&3\end{bmatrix},\quad f(t)=\begin{bmatrix}1\\-t\\t\end{bmatrix},\quad c=\begin{bmatrix}1\\0\\-1\end{bmatrix}.
   $$
   > [!IMPORTANT]
   > $A$ 的最小多项式 $m_A(\lambda)=(\lambda-2)^2$，令 $p_t(\lambda)=a_0(t)+a_1(t)\lambda$。由
   > $$
   > \begin{cases}a_0+2a_1=e^{2t}\\ a_1=te^{2t}\end{cases}
   > \Rightarrow \begin{cases}a_0=(1-2t)e^{2t}\\ a_1=te^{2t}\end{cases}
   > $$
   > 得 $e^{At}=a_0I+a_1A$，代入通解公式可求得完整解。

4. ==应用三：高阶常系数微分方程（状态方程）== 考察
   $$
   y^{(n)}+a_{n-1}y^{(n-1)}+\cdots+a_0y=u(t).
   $$
   令 $x_1=y,x_2=\dot{y},\dots,x_n=y^{(n-1)}$，则方程化为
   $$
   \dot{x}=Ax+f(t),\quad
   A=\begin{bmatrix}
   0&1&0&\cdots&0\\
   0&0&1&\cdots&0\\
   \vdots&\vdots&\vdots&\ddots&\vdots\\
   0&0&0&\cdots&1\\
   -a_0&-a_1&-a_2&\cdots&-a_{n-1}
   \end{bmatrix},\quad
   f(t)=\begin{bmatrix}0\\0\\\vdots\\0\\u(t)\end{bmatrix}.
   $$
   $A$ 称为系统矩阵，$x(t)$ 称为状态向量。这样即可用一阶方程组的公式求解。

5. ==应用四：系统的可控性== 对线性定常系统 $\dot{x}=Ax+Bu$，定义
   $$
   W(0,t_1)=\int_0^{t_1}e^{-At}BB^Te^{-A^Tt}dt,\quad
   C=[B,AB,\dots,A^{n-1}B].
   $$
   该系统完全可控的充要条件是 $W(0,t_1)$ 正定，或等价地 $\operatorname{rank}(C)=n$。

## 1.3 投影矩阵及其应用

### 1.3.1 投影与正交投影

1. ==定义3.1.1（投影）== 设 $W_1,W_2$ 是线性空间 $V$ 的两个子空间且 $V=W_1+W_2$，则对任意 $x\in V$ 可唯一分解为 $x=y+z$（$y\in W_1,z\in W_2$），称 $y$ 为 $x$ 沿 $W_2$ 在 $W_1$ 上的投影。

2. ==定义3.1.2（正交投影）== 若 $V=W_1\oplus W_2$（直和且正交），则称上述 $y$ 为 $x$ 在 $W_1$ 上的正交投影。

3. ==定理3.1.1（投影定理）== 设 $W$ 是酉空间 $V$ 的有限维子空间，则 $\forall x\in V$，$x$ 在 $W$ 上的正交投影存在且唯一。

4. ==例3.1.1== 设 $W_1=\operatorname{span}\{[1,0,1]^T\}$，$W_2=\{(x_1,x_2,x_3)^T\mid x_3=0\}$，求 $x=[1,1,1]^T$ 在 $W_2$ 上的投影。
   > [!IMPORTANT]
   > 取 $W_1$ 基 $\varepsilon_1=[1,0,1]^T$，$W_2$ 基 $\varepsilon_2=[1,0,0]^T,\varepsilon_3=[0,1,0]^T$，则 $\{\varepsilon_1,\varepsilon_2,\varepsilon_3\}$ 构成 $V$ 的基。设 $x=k_1\varepsilon_1+k_2\varepsilon_2+k_3\varepsilon_3$，解方程组得 $k_1=1,k_2=0,k_3=1$，故在 $W_2$ 上的投影为 $[0,1,0]^T$。

5. ==例3.1.2== 设 $W=\operatorname{span}(x_1,x_2)$，其中 $x_1=[2,5,-1]^T,x_2=[-2,1,1]^T$，求 $x=[1,2,3]^T$ 在 $W$ 和 $W^\perp$ 上的正交投影。
   > [!IMPORTANT]
   > 方法一：先求 $W^\perp$ 的基 $x_3=[7,0,14]^T$，再求 $x$ 在 $W^\perp$ 上的投影 $z=\frac{(x,x_3)}{(x_3,x_3)}x_3=[\frac{7}{5},0,\frac{14}{5}]^T$，则 $y=x-z=[-\frac{2}{5},2,\frac{1}{5}]^T$ 为 $W$ 上的投影。
   > 方法二：直接利用正交投影公式 $\operatorname{Proj}_W x=\frac{(x,x_1)}{\|x_1\|^2}x_1+\frac{(x,x_2)}{\|x_2\|^2}x_2$ 计算得相同结果。

6. ==命题3.1.1（正交投影的基表示）== 若 $x_1,\dots,x_n$ 是 $W$ 的一组正交基，则对任意 $y\in V$，
   $$
   \operatorname{Proj}_W y=\frac{(y,x_1)}{(x_1,x_1)}x_1+\cdots+\frac{(y,x_n)}{(x_n,x_n)}x_n.
   $$
   特别地，若为**标准正交基**，则 $\operatorname{Proj}_W y=(y,x_1)x_1+\cdots+(y,x_n)x_n$。

7. ==定理3.1.2（Gram-Schmidt正交化）== 有限维内积空间必存在标准正交基。构造方法：
   > [!IMPORTANT]
   > 设 $x_1,\dots,x_n$ 是一组基，令
   > $$
   > y_1=x_1,\quad
   > y_k=x_k-\sum_{i=1}^{k-1}\frac{(x_k,y_i)}{(y_i,y_i)}y_i,\quad
   > z_k=\frac{y_k}{\|y_k\|},
   > $$
   > 则 $z_1,\dots,z_n$ 即为标准正交基。进而有 $(x_1,\dots,x_n)=(z_1,\dots,z_n)A$，其中 $A$ 为正线上三角矩阵，对角线元素为 $\|y_i\|$。

### 1.3.2 最佳逼近

1. ==定义3.1.3（最佳逼近）== 设 $W$ 是线性空间 $V$ 的非空子集，$\alpha\in V$，若存在 $x\in W$ 满足
   $$
   \|\alpha-x\|\leq\|\alpha-y\|,\quad\forall y\in W,
   $$
   则称 $x$ 是 $\alpha$ 在 $W$ 中的最佳逼近向量。

2. ==定理3.1.3（最佳逼近定理）== 设 $W$ 是线性空间 $V$ 的线性子空间，则 $V$ 中任一向量 $x$ 在 $W$ 上都有唯一的最佳逼近，且该最佳逼近正是 $x$ 在 $W$ 上的正交投影。

3. ==定理3.1.4（最佳逼近的充要条件）== 设 $W$ 是酉空间 $V$ 的子空间，$x\in V$，则 $y\in W$ 是 $x$ 在 $W$ 上的最佳逼近当且仅当 $x-y\perp W$。
   > [!NOTE]
   > **必要性证明：** 设 $y$ 是最佳逼近但 $x-y\perp W$ 不成立，则存在 $z\in W$（$\|z\|=1$）使 $(x-y,z)=t\neq0$。令 $u=y+tz\in W$，则 $\|x-u\|^2=\|x-y\|^2-|t|^2<\|x-y\|^2$，矛盾。
   > **充分性证明：** 设 $x-y\perp W$，则对任意 $z\in W$，有$\|x-z\|^2=\|x-y\|^2+\|y-z\|^2\geq\|x-y\|^2$。

4. ==最佳逼近的代数表示== 设 $e_1,\dots,e_r$ 是 $W$ 的一组基，$y=k_1e_1+\cdots+k_re_r$ 是 $x$ 的最佳逼近。由 $x-y\perp e_i$ 得方程组
   $$
   \begin{bmatrix}
   (e_1,e_1) & \cdots & (e_r,e_1) \\
   \vdots & \ddots & \vdots \\
   (e_1,e_r) & \cdots & (e_r,e_r)
   \end{bmatrix}
   \begin{bmatrix}k_1\\ \vdots \\ k_r\end{bmatrix}
   =\begin{bmatrix}(x,e_1)\\ \vdots \\ (x,e_r)\end{bmatrix}.
   $$
   左侧矩阵为度量矩阵，必可逆，故 $k$ 唯一确定。若 $e_i$ 为标准正交基，度量矩阵为单位阵，则 $k_i=(x,e_i)$。

### 1.3.3 投影矩阵的应用

1. ==应用一：最小二乘问题（数据拟合）== 已知 $y$ 与 $x_1,\dots,x_n$ 满足线性关系 $y=k_1x_1+\cdots+k_nx_n$，通过 $m(\geq n)$ 次观测得数据，求系数 $k_i$ 使得残差平方和最小：
   $$
   \min_{l\in\mathbb{C}^n}\sum_{j=1}^m\Bigl|y^{(j)}-\sum_{i=1}^n l_i x_i^{(j)}\Bigr|^2.
   $$
   令 $a_i=(x_i^{(1)},\dots,x_i^{(m)})^T$，$b=(y^{(1)},\dots,y^{(m)})^T$，$A=(a_1,\dots,a_n)$，则问题转化为求 $b$ 在 $R(A)$ 上的最佳逼近，正规方程为 $A^HAk=A^Hb$。

2. ==例3.1.5（线性回归）== 已知 $y$ 与 $x$ 成线性关系 $y=ax+b$，观测数据为
   | $x$  | 3.6  | 3.7  | 3.8  | 3.9  | 4.0  | 4.1  | 4.2  |
   | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
   | $y$  | 1.00 | 0.90 | 0.90 | 0.81 | 0.60 | 0.56 | 0.35 |
   
   令 $x_1=a,x_2=b$，得
   $$
   A=\begin{bmatrix}3.6&1\\3.7&1\\3.8&1\\3.9&1\\4.0&1\\4.1&1\\4.2&1\end{bmatrix},\quad
   b=\begin{bmatrix}1\\0.9\\0.9\\0.81\\0.6\\0.56\\0.35\end{bmatrix},
   $$
   解 $A^TAx=A^Tb$ 得 $x=[-1.05,4.81]^T$，故 $y=-1.05x+4.81$。

3. ==应用二：函数的最佳逼近== 设 $f\in C[a,b]$，$\varphi_1,\dots,\varphi_n$ 为线性无关函数组，求$P(x)=\sum_{i=1}^n k_i\varphi_i(x)$ 使得 $\int_a^b(f(x)-P(x))^2dx$ 最小。
   
   > [!IMPORTANT]
   > 定义内积 $(f,g)=\int_a^b f(t)g(t)dt$，则问题转化为求 $f$ 在 $W=\operatorname{span}\{\varphi_1,\dots,\varphi_n\}$ 上的最佳逼近。$k$ 满足
   > $$
   > \begin{bmatrix}
   > (\varphi_1,\varphi_1) & \cdots & (\varphi_1,\varphi_n) \\
   > \vdots & \ddots & \vdots \\
   > (\varphi_n,\varphi_1) & \cdots & (\varphi_n,\varphi_n)
   > \end{bmatrix}
   > \begin{bmatrix}k_1\\ \vdots \\ k_n\end{bmatrix}
   > =\begin{bmatrix}(f,\varphi_1)\\ \vdots \\ (f,\varphi_n)\end{bmatrix}.
   > $$
   > 当 $\varphi_i$ 取多项式基或三角函数基时，即为常见的最小二乘曲线拟合或傅里叶逼近。


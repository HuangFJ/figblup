use na::DMatrix;

/// 个体间的加性遗传相关矩阵 $A$，其主对角线元素 $a_{ii}=1+F_i$，$F_i$ 是个体 $i$ 的近交系数。
/// $A$ 与加性遗传方差相乘 $A\sigma_u^2$ 是个体育种值间的协方差，如果 $u_i$ 是个体 $i$ 的育种值，则：
/// $$\begin{aligned}
///     var(u_i)
///     &= a_{ii}\sigma_u^2 \cr
///     &=(1+F_i)\sigma_u^2
/// \end{aligned}$$
///
/// 矩阵 $A$ 可以分解为：
/// $$
///     A = TDT^T
/// $$
///
/// $D$ 是对角矩阵，元素 $d_{ii}$ 是个体 $i$ 的孟德尔抽样的相关系数。
/// 如果个体的育种值 $u_i$ 与，亲本的育种值是 $u_s$ 和 $u_d$，则孟德尔抽样 $m_i=u_i-0.5(u_s+u_d)$ 。$D$ 的元素：
/// $$\begin{aligned}
///     d_{ii} &= var(m_i)/\sigma_u^2 \cr
/// 如果有双亲：d_{ii} &= 0.5-0.25(F_s+F_d) \cr
/// 如果只有单亲：d_{ii} &= 0.75-0.25(F_s) \cr
/// 双亲均无：d_{ii} &= 1
/// \end{aligned}$$
///
/// $T$ 是下三角矩阵，元素 $t_{ij}$ 是个体 $i$ 与个体 $j$ 的亲缘系数，并假设不存在近交，从一个世代到下一个世代追溯基因的流动，只需解释亲本到后代的关系。
/// $$\begin{aligned}
///     t_{ii} &= 1 \cr
/// 如果有双亲：t_{ij} &= 0.5(t_{sj}+t_{dj}) \cr
/// 如果只有单亲：t_{ij} &= 0.5(t_{sj}) \cr
/// 双亲均无：t_{ij} &= 0
/// \end{aligned}$$
///
/// 例子：
/// ```
/// let std_ped_vec = vec![
///     (1, 0, 0),
///     (2, 0, 0),
///     (3, 1, 2),
///     (4, 1, 0),
///     (5, 4, 3),
///     (6, 5, 2),
/// ];
/// let n = std_ped_vec.len();
/// let mut a = DMatrix::<f64>::zeros(n, n);
/// a_mat(&std_ped_vec, &mut a);
/// println!("{}", a.try_inverse().unwrap());
/// ```
pub fn a_mat(std_ped_vec: &Vec<(i64, i64, i64)>, mat: &mut DMatrix<f64>) {
    for ind in std_ped_vec {
        let i_idx = (ind.0 - 1) as usize;
        if ind.1 > 0 && ind.2 > 0 {
            //双亲都已知
            let s_idx = (ind.1 - 1) as usize;
            let d_idx = (ind.2 - 1) as usize;
            for j_idx in 0..i_idx {
                mat[(i_idx, j_idx)] = 0.5 * (mat[(j_idx, s_idx)] + mat[(j_idx, d_idx)]);
                mat[(j_idx, i_idx)] = mat[(i_idx, j_idx)];
            }
            mat[(i_idx, i_idx)] = 1.0 + 0.5 * mat[(s_idx, d_idx)];
        } else if ind.1 == 0 && ind.2 == 0 {
            //双亲均未知，且双亲不相关
            for j_idx in 0..i_idx {
                mat[(i_idx, j_idx)] = 0.0;
                mat[(j_idx, i_idx)] = mat[(i_idx, j_idx)];
            }
            mat[(i_idx, i_idx)] = 1.0;
        } else {
            //仅知一个亲本，并假设与配偶不相关
            let p = if ind.1 > 0 { ind.1 } else { ind.2 };
            let p_idx = (p - 1) as usize;
            for j_idx in 0..i_idx {
                mat[(i_idx, j_idx)] = 0.5 * mat[(j_idx, p_idx)];
                mat[(j_idx, i_idx)] = mat[(i_idx, j_idx)];
            }
            mat[(i_idx, i_idx)] = 1.0;
        }
    }
}

/// 迭代求解 $A^{-1}$ 矩阵。复杂度 O(n*n)
/// 例子：
/// ```
/// let std_ped_vec = vec![
///     (1, 0, 0),
///     (2, 0, 0),
///     (3, 1, 2),
///     (4, 1, 0),
///     (5, 4, 3),
///     (6, 5, 2),
/// ];
/// let n = std_ped_vec.len();
/// let mut a_inv = DMatrix::<f64>::zeros(n, n);
/// inv_a_mat(&std_ped_vec, &mut a_inv);
/// println!("{}", &a_inv);
/// ```
pub fn inv_a_mat(std_ped_vec: &Vec<(i64, i64, i64)>, mat: &mut DMatrix<f64>) {
    let n = std_ped_vec.len();
    // 如果令 $L=T\sqrt{D}$ ，则 $A=LL^T$
    let mut l_mat = DMatrix::<f64>::zeros(n, n);
    l_mat[(0, 0)] = 1.0;
    for &(ind, sire, dam) in std_ped_vec {
        let i_idx = (ind - 1) as usize;
        // $D^{-1}$ 的对角元素。$D$ 的对角元素 $d_{ii}=l_{ii}^2$
        let d_inv_ii = 1.0 / upd_l_mat(ind, sire, dam, &mut l_mat);

        mat[(i_idx, i_idx)] += d_inv_ii;

        if sire > 0 && dam > 0 {
            //双亲都已知
            let s_idx = (sire - 1) as usize;
            let d_idx = (dam - 1) as usize;

            //亲子加性
            let c_p_additive = -d_inv_ii / 2.0;
            mat[(i_idx, s_idx)] += c_p_additive;
            mat[(i_idx, d_idx)] += c_p_additive;
            mat[(s_idx, i_idx)] += c_p_additive;
            mat[(d_idx, i_idx)] += c_p_additive;

            //双亲加性
            let p_p_additive = d_inv_ii / 4.0;
            mat[(s_idx, s_idx)] += p_p_additive;
            mat[(d_idx, d_idx)] += p_p_additive;
            mat[(s_idx, d_idx)] += p_p_additive;
            mat[(d_idx, s_idx)] += p_p_additive;
        } else if sire > 0 || dam > 0 {
            //仅知一个亲本
            let p = if sire > 0 { sire } else { dam };
            //亲子加性
            let c_p_additive = -d_inv_ii / 2.0;
            mat[(i_idx, (p - 1) as usize)] += c_p_additive;
            mat[((p - 1) as usize, i_idx)] += c_p_additive;

            //双亲加性
            let p_p_addtioner = d_inv_ii / 4.0;
            mat[((p - 1) as usize, (p - 1) as usize)] += p_p_addtioner;
        }
    }
    // println!("{}", &l_mat);
}

/// 计算 $L$ 矩阵的第 $i$ 行元素，返回 $l_{ii}$ 的平方：
/// $$\begin{aligned}
///     l_{ii}
///     &=\sqrt{1.0-0.25(\sum_{m=1}^{s}l_{sm}^2+\sum_{m=1}^{d}l_{dm}^2)} \cr
///     l_{ij}
///     &=l_{ji} \cr
///     &=0.5(l_{sj}+l_{dj}) \hspace1em 其中 s 和 d 等于或大于 j
/// \end{aligned}$$
fn upd_l_mat(ind: i64, sire: i64, dam: i64, mat: &mut DMatrix<f64>) -> f64 {
    let i_idx = (ind - 1) as usize;
    for j_idx in 0..i_idx {
        let l_sj = if j_idx as i64 <= (sire - 1) {
            mat[((sire - 1) as usize, j_idx)]
        } else {
            0.0
        };
        let l_dj = if j_idx as i64 <= (dam - 1) {
            mat[((dam - 1) as usize, j_idx)]
        } else {
            0.0
        };
        mat[(i_idx, j_idx)] = 0.5 * (l_sj + l_dj);
    }

    let sum_s = (0..sire)
        .map(|m_idx| mat[((sire - 1) as usize, m_idx as usize)].powi(2))
        .reduce(|accum, item| accum + item)
        .unwrap_or(0.0);
    let sum_d = (0..dam)
        .map(|m_idx| mat[((dam - 1) as usize, m_idx as usize)].powi(2))
        .reduce(|accum, item| accum + item)
        .unwrap_or(0.0);

    let l_ii_sqr = 1.0 - 0.25 * (sum_s + sum_d);
    mat[(i_idx, i_idx)] = l_ii_sqr.sqrt();

    l_ii_sqr
}

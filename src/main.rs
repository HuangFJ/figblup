extern crate nalgebra as na;
// numerator relationship matrix A
//https://docs.rs/nalgebra/latest/nalgebra/base/struct.Matrix.html

mod gen_cov;
mod utils;

use gen_cov::*;
use na::*;
use utils::*;

/**
$$
\begin{bmatrix}
X^TX & X^TZ \cr
Z^TX & Z^TZ+\alpha A^{-1}
\end{bmatrix}
\begin{bmatrix}
\hat{\mu} \cr
\hat{u}
\end{bmatrix} =
\begin{bmatrix}
X^Ty \cr
Z^Ty
\end{bmatrix}
$$

$$
y_{ij}=\mu_i+u_j+e_ij
$$
$y_{ij}$ 是 $i$ 性别 $j$ 个体的 WWG；$\mu_i$ 是 $i$ 性别的固定效应；$u_j$ 是 $j$ 个体的随机效应；$e_{ij}$ 是随机误差效应。
假设：$\sigma_u^2 =20;\sigma_e^2 =40$ ，则：$\alpha =\sigma_e^2/\sigma_u^2=2$

<pre>
ID Sex BoarID  SowID WWG（kg）
4  公    1     Nil   4.5
5  母    3     2     2.9
6  母    1     2     3.9
7  公    4     5     3.5
8  公    3     6     5.0
</pre>
```
let std_ped_vec = vec![
    (1, 0, 0),
    (2, 0, 0),
    (3, 0, 0),
    (4, 1, 0),
    (5, 3, 2),
    (6, 1, 2),
    (7, 4, 5),
    (8, 3, 6),
];
let n = std_ped_vec.len();
let mut a_inv_mat = DMatrix::<f64>::zeros(n, n);
// 计算A的逆矩阵
inv_a_mat(&std_ped_vec, &mut a_inv_mat);
```
X 固定效应（性别）矩阵
<pre>
  公  母
4 1   0
5 0   1
6 0   1
7 1   0
8 1   0
</pre>
```
let x_mat = dmatrix![
    1.0,  0.0;
    0.0,  1.0;
    0.0,  1.0;
    1.0,  0.0;
    1.0,  0.0;
];
```
Z 是记录个体与所有个体的数据有无情况，不管有没有数据。1-3的个体是亲本，没有记录
<pre>
    1   2   3   4   5   6   7   8
4               1
5                   1
6                       1
7                           1
8                               1
</pre>
```
let z_mat = dmatrix![
    0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0;
    0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0;
    0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0;
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0;
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0;
];

let y_vec = dvector![4.5, 2.9, 3.9, 3.5, 5.0];

let x_t_x = &x_mat.transpose() * &x_mat;
let x_t_z = &x_mat.transpose() * &z_mat;
let z_t_x = x_t_z.transpose();
let z_t_z = &z_mat.transpose() * &z_mat;
let x_t_y = &x_mat.transpose() * &y_vec;
let z_t_y = &z_mat.transpose() * &y_vec;

let alpha = 2.0;
let r_mat = compose_sqr_mat(x_t_x, x_t_z, z_t_x, z_t_z + (a_inv_mat * alpha));

let y = compose_v_vec(x_t_y, z_t_y);

println!("{}", r_mat.try_inverse().unwrap() * y);
```
结果是：
<pre>
┌         ┐
│   4.358 │
│   3.404 │
│   0.098 │
│  -0.019 │
│  -0.041 │
│  -0.009 │
│  -0.186 │
│   0.177 │
│  -0.249 │
│   0.183 │
└         ┘
</pre>
*/
fn main() {
    let std_ped_vec = vec![
        (1, 0, 0),
        (2, 0, 0),
        (3, 0, 0),
        (4, 1, 0),
        (5, 3, 2),
        (6, 1, 2),
        (7, 4, 5),
        (8, 3, 6),
    ];
    let n = std_ped_vec.len();
    let mut a_inv_mat = DMatrix::<f64>::zeros(n, n);
    // 计算A的逆矩阵
    inv_a_mat(&std_ped_vec, &mut a_inv_mat);
    let x_mat = dmatrix![
        1.0,  0.0;
        0.0,  1.0;
        0.0,  1.0;
        1.0,  0.0;
        1.0,  0.0;
    ];
    let z_mat = dmatrix![
        0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0;
        0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0;
    ];
    let y_vec = dvector![4.5, 2.9, 3.9, 3.5, 5.0];

    let x_t_x = &x_mat.transpose() * &x_mat;
    let x_t_z = &x_mat.transpose() * &z_mat;
    let z_t_x = x_t_z.transpose();
    let z_t_z = &z_mat.transpose() * &z_mat;
    let x_t_y = &x_mat.transpose() * &y_vec;
    let z_t_y = &z_mat.transpose() * &y_vec;

    let alpha = 2.0;
    let r_mat = compose_sqr_mat(x_t_x, x_t_z, z_t_x, z_t_z + (a_inv_mat * alpha));

    let y = compose_v_vec(x_t_y, z_t_y);

    println!("{}", r_mat.try_inverse().unwrap() * y);
}

# Fast Computational Method

- operation_time:

    計算加法、乘法及三角函數(sin)所需要的時間。

- compute_matrix:

    N×N 矩陣 A 乘上一個 N×1 向量 x 的程式，透過 openmp 做平行計算。

- quick_median:

    從 quick sort 改寫出 quick median，假設一個長度為 N 的一個數列，找出中間的 N/2 大的數。

- fourier_transform:

    The discrete Fourier transform transforms a sequence of N complex numbers by
    
    ![formula].
        
    [formula]:https://wikimedia.org/api/rest_v1/media/math/render/svg/1af0a78dc50bbf118ab6bd4c4dcc3c4ff8502223
        
- fast_fourier_transform:

    以 bit-reverse 和 butterfly structure，對長度為 N 的陣列做DFT，N is radix-2, 3 and 5。
    
    利用 OpenMP 加速，平行處理每一層的所有 butterfly structure。

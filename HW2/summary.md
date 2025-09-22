# Fourier Modeling of Joint Angles and Stride Computation

## Overview
We model hip, knee, and ankle joint angles over the gait cycle using truncated Fourier series fit by linear least squares to all patient samples. This provides smooth, periodic representations suitable for downstream kinematic computations (stride distance and foot path).

## From FFT to a Truncated Fourier Model
Given a periodic signal $x(t)$ sampled over one gait cycle of length $N$, the real Fourier series with $K$ harmonics is
$$
\theta(t) \approx a_0 + \sum_{k=1}^{K} \big[a_k \cos(\omega_k t) + b_k \sin(\omega_k t)\big],\quad \omega_k = \tfrac{2\pi k}{N}.
$$
In our implementation, time is indexed by integer samples $x \in \{0,1,\dots,N-1\}$. The coefficients $a_0, a_k, b_k$ are estimated via linear least squares on a Fourier design matrix built from stacked patient data. An FFT of the per-time mean is computed only as an initialization/diagnostic; the final coefficients come from least-squares fitting to all raw points.

## Generic Discrete Model Used
Let $N$ be the number of samples per cycle. With $K$ harmonics, for sample index $x$ we have
$$
\theta(x) = a_0 + \sum_{k=1}^{K} \Big[a_k \cos\!\Big(\tfrac{2\pi k}{N} x\Big) + b_k \sin\!\Big(\tfrac{2\pi k}{N} x\Big)\Big].
$$
Coefficients are arranged as $[a_0, a_1, b_1, \dots, a_K, b_K]$, matching the summaries.

## Derived Joint Models (Degrees)
The fitted models below use the sample-index form with the dataset-specific $N$ and $K$.

### Hip (K = 4)
Coefficients: $a_0 = 14.52729527$; $(a_1, b_1) = (23.68382111, -3.454646754)$; $(a_2, b_2) = (-3.600895778, 0.7006881179)$; $(a_3, b_3) = (-0.2135900636, 1.734254375)$; $(a_4, b_4) = (-0.06123894065, 0.3105526648)$.

$$
\theta_{hip}(x) = 14.52729527 + \sum_{k=1}^{4} \Big[a_k^{(hip)} \cos\!\Big(\tfrac{2\pi k}{N_{hip}} x\Big) + b_k^{(hip)} \sin\!\Big(\tfrac{2\pi k}{N_{hip}} x\Big)\Big].
$$

See `hip_chart.png`: thin colored lines are individual patients, the bold black line is the sample-wise average, and the bold red line is the fitted Fourier equation evaluated from these coefficients.

![Hip Chart](hip_chart.png)

### Knee (K = 4)
Coefficients: $a_0 = 24.31908546$; $(a_1, b_1) = (-2.430398143, -18.24182655)$; $(a_2, b_2) = (-12.42330802, 10.88582658)$; $(a_3, b_3) = (-0.1963700052, 3.994124334)$; $(a_4, b_4) = (-0.5319396146, 0.4237436464)$.

$$
\theta_{knee}(x) = 24.31908546 + \sum_{k=1}^{4} \Big[a_k^{(knee)} \cos\!\Big(\tfrac{2\pi k}{N_{knee}} x\Big) + b_k^{(knee)} \sin\!\Big(\tfrac{2\pi k}{N_{knee}} x\Big)\Big].
$$

See `knee_chart.png`: thin colored lines are individual patients, the bold black line is the sample-wise average, and the bold red line is the fitted Fourier equation evaluated from these coefficients.

![Knee Chart](Knee_chart.png)

### Ankle (K = 5)
Coefficients: $a_0 = -5.67336993$; $(a_1, b_1) = (0.864044135, 6.418378396)$; $(a_2, b_2) = (-2.642263687, -7.277654774)$; $(a_3, b_3) = (-2.360690859, 2.640891244)$; $(a_4, b_4) = (1.712862189, -1.351434239)$; $(a_5, b_5) = (-0.5338225615, -0.1974550647)$.

$$
\theta_{ankle}(x) = -5.67336993 + \sum_{k=1}^{5} \Big[a_k^{(ankle)} \cos\!\Big(\tfrac{2\pi k}{N_{ankle}} x\Big) + b_k^{(ankle)} \sin\!\Big(\tfrac{2\pi k}{N_{ankle}} x\Big)\Big].
$$

See `ankle_chart.png`: thin colored lines are individual patients, the bold black line is the sample-wise average (after the +15° offset applied per the assignment), and the bold red line is the fitted Fourier equation evaluated from these coefficients.

![Ankle Chart](ankle_chart.png)

Note: The ankle series was fit after adding +15° to patient ankle signals per the assignment; the fitted function reflects that offset.

## Gait Path (Foot Position) from Joint Angles
Let segment lengths (in cm) be $L_{top} = 46$, $L_{bottom} = 45$, and $A = 9.4$. Define
$$
\phi_{hip}(x) = (\theta_{hip}(x) - 90)\tfrac{\pi}{180},\quad
\phi_{knee}(x) = \theta_{knee}(x)\tfrac{\pi}{180},\quad
\phi_{ankle}(x) = \theta_{ankle}(x)\tfrac{\pi}{180}.
$$
The foot position relative to the hip is
$$
\begin{aligned}
X(x) &= L_{top}\cos\phi_{hip}(x) + L_{bottom}\cos\big(\phi_{hip}(x) - \phi_{knee}(x)\big) + A\cos\big(\phi_{hip}(x) - \phi_{knee}(x) + \phi_{ankle}(x)\big),\\
Y(x) &= L_{top}\sin\phi_{hip}(x) + L_{bottom}\sin\big(\phi_{hip}(x) - \phi_{knee}(x)\big) + A\sin\big(\phi_{hip}(x) - \phi_{knee}(x) + \phi_{ankle}(x)\big).
\end{aligned}
$$
These are evaluated on a common timeline by reparameterizing each joint’s Fourier series to a shared length $N_c$ and sampling $\theta_{\cdot}(x)$ at $x = 0,1,\dots,N_c-1$. See `stride_path.png` for the parametric $(X,Y)$ trajectory over a cycle.

![Stride Path](stride_path.png)

## Stride Length Approximation
The stride distance along the horizontal direction over one cycle is approximated as the span of $X(x)$:
$$
\text{StrideLength} \approx \max_x X(x) - \min_x X(x).
$$
This matches the computation that yielded the reported values in `stride_summary.txt`:
- $X_{min} = -44.290294$ cm, $X_{max} = 48.773133$ cm → Full stride $\approx 93.063427$ cm.

See `stride_chart.png` for the $X(x)$ time series used in this calculation (blue line).

![Stride Chart](stride_chart.png)

A more refined definition could integrate arc length along the $(X, Y)$ foot path or measure heel-strike-to-heel-strike displacement with stance/swing segmentation; that would require additional modeling assumptions and potentially modified linkage geometry or event detection, hence “new equations based on an altered gait path.”

## Velocity and Angular Frequency
If we scale the temporal frequency (cadence), letting $x\mapsto x$ but replacing $\omega_k = 2\pi k/N$ by $\tilde{\omega}_k = \alpha\,\omega_k$ with $\alpha>1$, we obtain a faster oscillation:
$$
\tilde{\theta}(x) = a_0 + \sum_{k=1}^{K} \Big[a_k \cos(\alpha\,\omega_k x) + b_k \sin(\alpha\,\omega_k x)\Big].
$$
The instantaneous foot velocity components follow by differentiation:
$$
\dot{X} = -L_{top}\sin\phi_{hip}\,\dot{\phi}_{hip} - L_{bottom}\sin(\phi_{hip}-\phi_{knee})\,(\dot{\phi}_{hip}-\dot{\phi}_{knee}) - A\sin(\phi_{hip}-\phi_{knee}+\phi_{ankle})\,(\dot{\phi}_{hip}-\dot{\phi}_{knee}+\dot{\phi}_{ankle}),
$$
with analogous $\dot{Y}$. Since $\dot{\phi}$ scales linearly with frequency, increasing $\alpha$ increases linear speed proportionally (for the same angular amplitude). Thus, gait velocity can be adjusted directly by increasing the angular frequency without refitting amplitudes, while changing stride length meaningfully would require altering the joint angle waveforms or the linkage geometry (the “new equations”).

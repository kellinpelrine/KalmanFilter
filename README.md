# KalmanFilter

Author: Kellin Pelrine

Code for a Kalman filter, Kalman smoother, and missing data interpolation.

Kalman_filter_smoothing.m illustrates a Kalman filter and Kalman smoother. It produces plots like the following two:

![image](https://github.com/kellinpelrine/KalmanFilter/blob/master/Filter%20Graph.jpg)

Here blue is the underlying state (not observed), red is the observations, and green is the filter predictions. 

![image](https://github.com/kellinpelrine/KalmanFilter/blob/master/Smoother%20Graph.jpg)

Blue is the state, red is the observations, black is the smoothed estimate (best prediction in hindsight).

Kalman_missing_data.m illustrates the use of a Kalman filter and Kalman smoother to impute missing data. It produces a plot like this:

![image](https://github.com/kellinpelrine/KalmanFilter/blob/master/Missing%20Data%20Graph.jpg)

Here blue is the state, red is the state where data is missing, green is the filter extrapolation for the missing portion, and black is the smoother interpolation for the missing portion.

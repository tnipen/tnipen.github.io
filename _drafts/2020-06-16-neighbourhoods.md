These methods can be computationally expensive. A naiive method has  computational complexity of O(GR^2), where G is the number of gridpoints and R is the full-width of the neighbourhood. Gridpp uses several tricks and approximations to improve the computational speed. To improve speed, the main strategy is to find ways that reuse calculations that have been done for one neighbourhood when processing another neighbourhood.

To find the neighbourhood **mean**, the domain is accumulated from the lower left corner to the upper right. This algorithm is O(G) and is independent of the neighbourhood size. A neighbourhood mean can be computed by the values at the 4 corners of the neighbourhood:
```python
v = (upper_left + lower_left - upper_right - lower_right) / R / R
```

To find the neighbourhood **standard deviation** or **variance**, the following relationship is used:
```
variance = E[x^2] - E[x]^2
```
That is, the variance is efficiently compute by calculating two types of neighbourhood means and taking the difference.

To find the neighbourhood **minimum** or **maximum**, the domain is split into Nx1 slices. A neighbourhood is then represented by R of these Rx1 slices. The minimum (maximum) is computed for each slice, and then a minimum (maximum) of R of the Rx1 slices is computed. This saves time because the Nx1 minimums (maximums) are reused R times. This algorithm has time complexity O(GRlog(R)).

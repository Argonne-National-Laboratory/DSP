bool isapprox(double lb, double mid, double ub) {
    if (lb > mid + 1.e-6) {
        return false;
    }
    if (ub < mid - 1.e-6) {
        return false;
    }
    return true;
}

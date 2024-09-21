#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>

const double EPS = 1e-9;
const int INF = 2; // it doesn't actually have to be infinity or a big number

int gauss (std::vector<std::vector<double>> a, std::vector<double> &ans) {
    int n = (int) a.size();
    int m = (int) a[0].size() - 1;

    std::vector<int> where (m, -1);
    for (int col=0, row=0; col<m && row<n; ++col) {
        int sel = row;
        for (int i=row; i<n; ++i)
            if (std::abs(a[i][col]) > std::abs(a[sel][col]))
                sel = i;
        if (std::abs(a[sel][col]) < EPS)
            continue;
        for (int i=col; i<=m; ++i)
            std::swap(a[sel][i], a[row][i]);
        where[col] = row;

        for (int i=0; i<n; ++i)
            if (i != row) {
                double c = a[i][col] / a[row][col];
                for (int j=col; j<=m; ++j)
                    a[i][j] -= a[row][j] * c;
            }
        ++row;
    }

    ans.assign(m, 0);
    for (int i=0; i<m; ++i)
        if (where[i] != -1)
            ans[i] = a[where[i]][m] / a[where[i]][i];
    for (int i=0; i<n; ++i) {
        double sum = 0;
        for (int j=0; j<m; ++j)
            sum += ans[j] * a[i][j];
        if (std::abs(sum - a[i][m]) > EPS)
            return 0; // No solution
    }

    for (int i=0; i<m; ++i)
        if (where[i] == -1)
            return INF; // Infinite solutions
    return 1; // Unique solution
}

int main() {
    // Example augmented matrix for a system of 2 equations and 2 variables:
    // x +  y = 2
    // 2x + 3y = 5
    std::vector<std::vector<double>> a = {
        {1, 1, 2},
        {2, 3, 5}
    };
    std::vector<double> ans;

    int result = gauss(a, ans);

    if (result == 0) {
        std::cout << "No solution\n";
    } else if (result == INF) {
        std::cout << "Infinite solutions\n";
    } else {
        std::cout << "Unique solution:\n";
        for (double x : ans) {
            std::cout << x << " ";
        }
        std::cout << "\n";
    }

    return 0;
}

#include <cstdint>
#include <iostream>
#include <cmath>

using namespace std;

constexpr size_t get_split(size_t m, size_t aggr) {
	size_t rounded_down = aggr * ((m / 2) / aggr);
	size_t rounded_up = aggr * ((m / 2 + aggr - 1) / aggr);
	return rounded_up < (m - rounded_down) ? rounded_up : rounded_down;
}

constexpr size_t get_split_fast(size_t m, size_t aggr) {
	return aggr * ((m / 2 + (aggr) / 2) / aggr);
}

int main() {
	for (size_t aggr = 1; aggr < 1500; ++aggr) {
		for (size_t m = aggr + 1; m < 3000; ++m) {
			int s = (int) get_split(m, aggr);
			int f = (int) get_split_fast(m, aggr);
			if (abs(s - (int)m/2) < abs(f - (int)m/2)) {
				cout << m << "," << aggr << ": " << s << "," << f << std::endl;
			}
		}
	}
}
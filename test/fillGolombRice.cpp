#include <array>
#include <iostream>

#include "src/function/AbstractParallelRecSplit.hpp"

using namespace sux::function;

template<int LEAF_SIZE>
void print_golomb_rice() {
	array<uint32_t, MAX_BUCKET_SIZE> result = fill_golomb_rice<LEAF_SIZE>();
	std::cout << "{";
	for (int i = 0; i < MAX_BUCKET_SIZE - 1; ++i) {
		if (i > 0 && i % 10 == 0)
			std::cout << "\n ";
		printf("0x%08X,", result[i]);
		if (i % 10 != 9)
			std::cout << " ";
	}
	printf("0x%08X\n}", result[MAX_BUCKET_SIZE - 1]);
	if (LEAF_SIZE < 24)
		std::cout << ",";
	std::cout << "\n";
}

template<int FROM, int TO>
void print_all_golomb_rice() {
	print_golomb_rice<FROM>();
	if constexpr (FROM < TO)
		print_all_golomb_rice<FROM + 1, TO>();
}

int main() {
	constexpr int fromLeafSize = 2;
	constexpr int toLeafSize = 24;
	std::cout << "static constexpr std::array<std::array<uint32_t, " << MAX_BUCKET_SIZE
		<< ">, " << (toLeafSize - fromLeafSize + 1) << "> golomb_rice_memo{{\n";
	print_all_golomb_rice<fromLeafSize, toLeafSize>();
	std::cout << "}};" << std::endl;
}
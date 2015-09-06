#include 


void histo_effsigma()
 std::ifstream file("program.txt");
    for(std::string word; file >> word; )
        std::cout << word << '\n';

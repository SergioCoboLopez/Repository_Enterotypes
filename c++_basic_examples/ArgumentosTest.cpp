// Program: MyArgs
#include <iostream>
#include <string>
#include <sstream> // for std::stringstream
#include <cstdlib> // for exit()

int main(int argc, char *argv[])
{
  std::cout << "There are " << argc << " arguments:\n";

  for (int count=0; count < argc; ++count)
    std::cout << count << " " << argv[count] << '\n';


  if (argc <= 1)
    {

      if (argv[0])
	std::cout << "Usage: " << argv[0] << " <number>" << '\n';
      else
	std::cout << "Usage: <program name> <number>" << '\n';

      exit(1);
    }

  std::stringstream convert0(argv[1]), convert1(argv[2]), convert2(argv[3]);

  int myint0; int myint1; int myint2;
  if (!(convert0 >> myint0)) // do the conversion
    myint0 = 0; // if conversion fails, set myint to a default value

  if (!(convert1 >> myint1)) 
    myint1 = 0;

    if (!(convert2 >> myint2)) 
    myint2 = 0; 

  std::cout << "Got integers: " << myint0 << '\n';

  
  return 0;
}

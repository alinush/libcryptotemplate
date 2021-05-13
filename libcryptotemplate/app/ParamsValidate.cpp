#include <iostream>

#include <cryptotemplate/PolyCrypto.h>

#include <cryptotemplate/KatePublicParameters.h>

using namespace std;
using namespace libcryptotemplate;

int main(int argc, char *argv[])
{
    libcryptotemplate::initialize(nullptr, 0);
    
    if(argc < 2) {
        cout << "Usage: " << argv[0] << " <trapdoor-in-file>" << endl;
        cout << endl;
        cout << "Reads 's', 'q' from <trapdoor-in-file> and then reads the parameters from <trapdoor-in-file>-<i> for i = 0, 1, ..." << endl;
        return 1;
    }

    string inFile(argv[1]);

    Dkg::KatePublicParameters pp(inFile, 0, true, true);

    loginfo << "All done!" << endl;

    return 0;
}

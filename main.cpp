/*
    $Id: main.cpp,v 1.4 2004/05/21 05:50:13 taku-ku Exp $;
 
   Copyright (C) 2004 Taku Kudo, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.
  
   This program is free software; you can redistribute it and/or modify
     it under the terendpoint of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
    
   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA
*/
#include <unistd.h>
#include <fstream>
#include "fanta.h"

#define OPT " [-m sigma] [-d] [-e] [-w] "

void usage (void)
{
    std::cout << "fanta implementation by Yilin Chen" << std::endl;
    std::cout << std::endl;
    std::cout << "options" << std::endl;
    std::cout << "  -h, show this usage help" << std::endl;
    
    std::cout << "methods" << std::endl;
    std::cout << "  -P, pruning employed, 0=none, 1=+check, 2=+prune" << std::endl;
    std::cout << "  -S, share computation" << std::endl;
    std::cout << "  -D, support definition, 0=exact, 1=expected, 2=probabilistic (default: 0)"<<std::endl;
    
    std::cout << "parameters" << std::endl;
    std::cout << "  -m, max patten length (default: 0, infinite)" << std::endl;
    std::cout << "  -n, max idential labeled node (default: 0, infinite)" << std::endl;
    std::cout << "  -e, epsilon, set the absolute error (default: 0.1)" << std::endl;
    std::cout << "  -d, delta, set the confidence (default: 0.1)" << std::endl;
    std::cout << "  -s, sigma, set the minimum support (default: 1.0)" << std::endl;
    std::cout << "  -p, phi, set the threshold for probabilistic (default: 1)" << std::endl;
    std::cout << "  -i, input, set the path of input (default: 'in')" << std::endl;
    std::cout << "  -o, output, set the path of output (default: 'out')" << std::endl;
    std::cout << std::endl;
}

int main (int argc, char **argv)
{
	unsigned int maxpat = 0xffffffff;
	unsigned int maxnode = 0xffffffff;
    unsigned int definition = 0;
    unsigned int phi = 1;
    unsigned int pruning = 0;
    double epsilon = 0.1;
    double delta = 0.1;
    double sigma = 1;
    bool share = false;
    
    char infile[20] = "in";
    char outfile[20] = "out";
	int opt;
	while ((opt = getopt(argc, argv, "P:SD:m:n:e:d:s:p:i:o:h")) != -1) {
        switch(opt) {
            case 'P':
                pruning = atoi(optarg);
                break;
            case 'S':
                share = true;
                break;
            case 'D':
                definition = atoi(optarg);
                break;
            case 'm':
                maxpat = atoi (optarg);
                break;
            case 'n':
                maxnode = atoi (optarg);
                break;
            case 'e':
                epsilon = atof(optarg);
                break;
            case 'd':
                delta = atof(optarg);
                break;
            case 's':
                sigma = atof(optarg);
                break;
            case 'p':
                phi = atoi(optarg);
                break;
            case 'i':
                sprintf(infile,"%s",optarg);
                break;
            case 'o':
                sprintf(outfile,"%s",optarg);
                break;
            case 'h':
            default:
                usage ();
                return -1;
        }
	}
    std::ifstream in(infile);
    std::ofstream out(outfile);
	FANTA::Fanta fanta;
    fanta.run(pruning, share, definition, maxpat, maxnode, epsilon, delta, sigma, phi, in, out);
}

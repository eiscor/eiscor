#!/bin/bash

cd ../..                                                                                                        
make install                                                                                                    
cd example/double 
rm example_d_symtrid_qr_race
make example_d_symtrid_qr_race
./example_d_symtrid_qr_race

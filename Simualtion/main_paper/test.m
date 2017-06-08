clear;clc

Hcap = ( 1/sqrt(2) ) * ( randn(2,2) + 1i*randn(2,2) );

Htil = [ ( 1/sqrt(2) ) * ( randn(2,2,10) + 1i*randn(2,2,10) ) ]

H = Hcap + Htil;
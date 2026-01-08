#! /bin/awk -f
BEGIN{
avm=0
asm=0
m_4=0
binder=0
d=0
c_r=0

}

{
avm+=$2

asm+=$3

m_4+=$4

binder+=$5
d+=$5^2

c_r+=$6
}

END{
binder=binder/NR
d=d/NR
delta=(d-binder^2)/(NR-1)

printf "%5.3f %18.12f %18.12f %18.12f %18.12f %18.12f %18.12f\n",$1,avm/NR,asm/NR,m_4/NR,binder,delta^0.5,c_r/NR
} 


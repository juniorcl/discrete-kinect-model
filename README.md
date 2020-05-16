# Modelos Cinéticos Discretos
Modelos utilizados para estudo de crescimento de superfícies. Para estudo de futuros alunos ou pessoas que estejam interessados em estudos relacionados.

## Sobre
Este é um repositório onde se encontram os programas utilizados durante a realização do trabalho de conclusão de curso (TCC) do Instituto Federal de ciência, educação Tecnologia Fluminense (IFF). O TCC foi apresentado como forma de conclusão do curso de Ciências da Natureza - Licenciatura em Física. O trabalho original pode ser acessado por meio deste [link](http://bd.centro.iff.edu.br/jspui/handle/123456789/1360) ou pelo portal do [IFF](http://portal1.iff.edu.br/nossos-campi/campos-centro).

## Modelos Abordados
- ** Deposição Aleatória (DA) **
Esse é o modelo mais simples. O sítio para deposição é escolhido aleatoriamente para deposição. Apesar de ser o modelo mais simples, este pode servir de base para os demais.

- __Deposição Aleatória com Relaxação (DAR)__
DAR é uma variação do modelo DA. Neste modelo a escolha do sítimo ocorre do mesmo jeio como o DA. Porém a partícula pode descrever uma trajetória pelo sítio até encontrar um sítio de menor altura.

- ** Deposição Balística (DB)
A deposição ocorre de forma semelhante aos modelos anteriores. Mas se a partícula que for depositação se encontrar com outra lateralmente, esta por sua vez se tornará parte do filme.

- ** Deposição Sólido sobre Sólido Restrito (RSOS) **
Neste modelo a deposição só ocorre se a atura de seu vizinho for de uma determinada altura. Caso contrário a partícula não se deposita.

- ** Wolf-Villan (WV) e Das Sarma e Tamboronea (DT) **
O modelo WV tem o objetivo de maximizar o número de coordenações, então a partícula se agregará preferencialmente em uma posição que seja, em geral, mais baixa do que seus vizinhos laterais. Enquanto que DT não há a preferência de maximização do número de coordenações das partículas como no modelo anterior, neste modelo apenas almeja-se aumentar o número de coordenações.

## Como usar?
Os modelos estão escritos na linguagem de programação Fortran. Durante o trabalho os compiladores utilizados foram ifort da intel e gfortran e o sistema operacional GNU/Linux. De forma simples, só é preciso baixar os arquivos .f90 e compilar como exemplificado abaixo.

> $ gfortran DA.f90

Posteriormente será preciso executar o código com o comando:

> $ ./a.out

## Licença
Este projeto é licensiado sobre a licença MIT - veja o arquivo [LICENSE.md](LICENSE.md) para mais detalhes.




























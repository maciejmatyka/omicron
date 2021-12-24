# omicron 

This simple model shows the difference in transmission rates between alpha, delta and omicron. Please do not take this as a serious scientific proof of any kind. I used simple model of bouncing circles with flags meaning "healthy" and "ill" where "healthy" may become "ill" with some probability P if they meet together. Now, depends on the variant of the virus the P changes.

Values of P that I took are relative, not definite. What I mean is that the goal was to see the difference between virus variants, not definite answer on how fast it spreads or anything like that. I calculated values of P using three sources which gave me
delta = 125% of alpha 
omicron = 4.2 x delta
Sources:
1) https://fortune.com/2021/12/08/omicron-covid-variant-data-more-transmissible-than-delta-new-study/
2) https://www.weforum.org/agenda/2021/11/what-makes-the-delta-variant-different-covid-19/
3) https://fortune.com/2021/12/08/omicron-covid-variant-data-more-transmissible-than-delta-new-study/

# About the code 

Source code for the simulation: https://github.com/maciejmatyka/

The code is dirty and "working" meaning that I didn't clean it for you. You may find it interesting to dig a bit and find vaccination, birth, death, lockdown and other mechanisms related to pandemia and virus spreading. Have fun and
keep in mind - this is a toy model only (but somehow informative on the other hand).

# Comment from author
The code is dirty and "working" meaning that I didn't clean it for you. 
You may find it interesting to dig a bit and find vaccination, birth, death, 
lockdown and other mechanisms related to pandemia and virus spreading. Have fun and
keep in mind - this is a toy model only (but somehow informative on the other hand).

# How it works?
The model distributes N agents randomly. They interact (collide) with each other. They move on straight trajectories between collisions. When they collide - they may exchange virus with some probability P, so ILL can affect healthy one. Innitially all are healthy and we infect one person at timestep n=100.


Music by F. Chopin, Fantaisie - Impromptu, Op. 66, perf. Frank Levy (musopen.org)
Icons are by Viewbaa Ratipat free from ttps://xnimrodx.gumroad.com/l/lwfIu
============================================
#covid19 #omicron #simulation

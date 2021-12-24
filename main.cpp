/*
*Toy model* of virus spreading simulation 
by Maciej Matyka 2021-12-24
http://panoramx.ift.uni.wroc.pl/~maq/eng/

Title: How omicron overtakes delta and alpha (toy model)

Model:
This simple model shows the difference in transmission rates 
between alpha, delta and omicron. Please do not take this as a 
serious scientific proof of any kind. I used simple model of 
bouncing circles with flags meaning "healthy" and "ill" where 
"healthy" may become "ill" with some probability P if they meet 
together. Now, depends on the variant of the virus the P changes.
Values of P that I took are relative, not definite. What I mean 
is that the goal was to see the difference between virus 
variants, not definite answer on how fast it spreads or anything 
like that. I calculated values of P using three sources which 
gave me
delta = 125% of alpha
omicron = 4.2 x delta
Sources:
1) https://fortune.com/2021/12/08/omicron-covid-variant-data-more-transmissible-than-delta-new-study/
2) https://www.weforum.org/agenda/2021/11/what-makes-the-delta-variant-different-covid-19/
3) https://fortune.com/2021/12/08/omicron-covid-variant-data-more-transmissible-than-delta-new-study/


The code was used for the video: 

Comment from author
The code is dirty and "working" meaning that I didn't clean it for you. 
You may find it interesting to dig a bit and find vaccination, birth, death, 
lockdown and other mechanisms related to pandemia and virus spreading. Have fun and
keep in mind - this is a toy model only (but somehow informative on the other hand).

How it works?
The model distributes N agents randomly. They interact (collide) with each other.
They move on straight trajectories between collisions. When they collide - they may
exchange virus with some probability P, so ILL can affect healthy one. 
Innitially all are healthy and we infect one person at timestep n=100.

*/

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
using namespace std;
#include <GL/glew.h>
#include <GL/glut.h>
#include <random> 
#include "CVector3.h"

#define WIDTH W
#define HEIGHT H
   
const int W = 1200;
const int H = 1200;

enum state // the code is more general, in the video I used only the HEALTHY/ILL mode
{
	HEALTHY,
	ILL,
	CONVALESCENT,
	VACCINATED,
	PASSEDAWAY
};

// random numbers
std::mt19937 mt(1000);//time(NULL));
std::uniform_real_distribution<float> ran01(0, 1);
std::uniform_real_distribution<float> ran0505(-0.5, 0.5);

// some parameters
float dt = 1.0;
int N = 5500;//100;//1200;
float R = 0.004;
#define ALPHAP 0.2					// https://fortune.com/2021/12/08/omicron-covid-variant-data-more-transmissible-than-delta-new-study/
#define DELTAP (ALPHAP+1.25*ALPHAP)	// https://www.weforum.org/agenda/2021/11/what-makes-the-delta-variant-different-covid-19/
#define OMICRONP (DELTAP*4.2)		// https://fortune.com/2021/12/08/omicron-covid-variant-data-more-transmissible-than-delta-new-study/
float PSPREAD = OMICRONP;//0.5;//5;//0.1;//0.2;
float PPASSEDAWAY = 0.12;//5;//0.1;//0.2;		// albo ozdrowieniec, albo umiera
float PBIRTH = 0.3;//urodziny
float PLOCKDOWNFLIGHT = 0.01;			// ruch w trakcie lockdownu
int step = 0;
#define CIRCLER 0.48	


float PROCSZCZEP = 0;//0.2;

int lockdown = 0;

const float timetoconvalescent = 125.0;
const float convalescenttohealthytime = 125.0;
class ball
{
public:
	ball(CVector3 _r, CVector3 _v, float _R, int _s): r(_r), v(_v), s(_s),R(_R),changet0(0) {}
	CVector3 r;
	CVector3 v;
	float R;
	int s;		// state
	float changet0;	// time
};
vector<ball> kulki;

void init(void)
{
	for(int i=0; i<N; i++)
	{
		const float DV = 0.001;
		float vx = ran0505(mt)*DV;
		float vy = ran0505(mt)*DV;
		float _R = R;//+0.1*ran01(mt)*R;	// slight variation
		//float _R = ran01(mt)*R;	// slight variation
		float x,y;
		int colliding;
		//cout << i << endl;
		float d;
		do // find empty spot
		{	
			colliding=0;
			x = ran01(mt);
			y = ran01(mt);
			for(int j=0; j<kulki.size(); j++)
				if((CVector3(x,y,0) - kulki[j].r).getLength() < _R )
					colliding = 1;
			d = (CVector3(x,y,0)-CVector3(0.5,0.5,0)).getLength();
		} while(colliding==1 || d>CIRCLER);

		kulki.push_back(ball(CVector3(x,y,0),CVector3(vx, vy,0),_R,HEALTHY));

		if(ran01(mt)<PROCSZCZEP)	// zaszczepiony
			kulki[kulki.size()-1].s = VACCINATED;
	}
}

void redisplay(int w, int h) 
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluOrtho2D(0,1,0,1);
	glMatrixMode(GL_MODELVIEW);
}
//KACPER

void draw()
{
	glLoadIdentity();
	gluLookAt(0, 0, 1, 0, 0, -1, 0, 1, 0);
	glEnable( GL_POINT_SMOOTH );
	glColor4f(0.7,0.7,0.7,0.9);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPointSize(10);
	float r,g,b;
	for(int i=0; i<kulki.size(); i++)
	{
			if(kulki[i].s==HEALTHY)				glColor4f(245/255.0,255/255.0,250/255.0,255/255.0);
			else if(kulki[i].s==CONVALESCENT)	glColor4f(120/255.0,190/255.0,130/255.0,255/255.0);
			else if(kulki[i].s==ILL)			glColor4f(245/255.0,120/255.0,140/255.0,255/255.0);
			else if(kulki[i].s==VACCINATED)		glColor4f(245/255.0,255/255.0,250/255.0,255/255.0);
		//color(0.5,0.5,0.5,0.5,0.5,0.5,0.5,1.0,1.0,0.5,0.80,0.90,0.30, r,g,b);
		//glColor4f(r,g,b,0.97);

//		if(kulki[i].s==ILL)
		glBegin(GL_POINTS);
			glVertex2f(kulki[i].r.x , kulki[i].r.y);
			//glVertex2f(0.5f,0.5f);
		glEnd();

	}

}

void render(void) 
{
	glClearColor(0.04,0.02,0.021,1);
	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
	
	draw();

	glutSwapBuffers();
	glutPostRedisplay();

}

void move()
{
	for(int i=0; i<kulki.size(); i++)
	if(!(kulki[i].s==PASSEDAWAY))
	{
		//if(!(i%10))
		//	continue;
		kulki[i].r = kulki[i].r + dt * kulki[i].v;
		/*if(kulki[i].r.x < 0) kulki[i].r.x = 1+kulki[i].r.x;
		if(kulki[i].r.y < 0) kulki[i].r.y = 1+kulki[i].r.y;
		if(kulki[i].r.x > 1) kulki[i].r.x = kulki[i].r.x-1;
		if(kulki[i].r.y > 1) kulki[i].r.y = kulki[i].r.y-1;*/
		//int CIRCLER = W;

		float d = (kulki[i].r-CVector3(0.5,0.5,0)).getLength();
		if( d>CIRCLER )
		{
			//kulki[i].r = kulki[i].r - dt * kulki[i].v;
			CVector3 de = kulki[i].r-CVector3(0.5,0.5,0);
			CVector3 n = de;
			n.doNormalize();
			CVector3 Fk =  -1.0f * n*(fabs(de.getLength()-CIRCLER));		// Hooke
			kulki[i].v = kulki[i].v + Fk*dt;
			/*// reflect r
			CVector3 n = CVector3(0.5,0.5,0)-kulki[i].r;
			n.doNormalize();
			CVector3 vn = n*kulki[i].v;
			CVector3 vt = kulki[i].v - vn;
			kulki[i].v = kulki[i].v - 2.0*vn;*/
			
		}
	}


}

//https://introcs.cs.princeton.edu/java/assignments/collisions.html
void collision()
{
// collision
	double m1;
	double m2;
	m1 = m2 = 1;
	// collisions
	// https://en.wikipedia.org/wiki/Elastic_collision
	for (int i = 0; i < kulki.size(); i++)
	if(kulki[i].s!=PASSEDAWAY)
	{
		for (int j = i+1; j < kulki.size(); j++)
		if(kulki[j].s!=PASSEDAWAY)
		{
			double sum1,sum2;
			sum1 = kulki[i].v.getLength() + kulki[j].v.getLength();
			
			CVector3 r1 = kulki[i].r;
			CVector3 r2 = kulki[j].r;

			// o-t--o / separate first (to time of contact)
			CVector3 t = kulki[i].r-kulki[j].r;
			double d = t.getLength();

			if(d < (kulki[i].R+kulki[j].R))
			{	
				// virus

				if(ran01(mt)<PSPREAD)
				{
					if(kulki[i].s==ILL && kulki[j].s==HEALTHY) 
					{
						kulki[j].s=ILL;
						kulki[j].changet0=step*dt;
					}
					if(kulki[j].s==ILL && kulki[i].s==HEALTHY)
					{
						kulki[i].s=ILL;
						kulki[i].changet0=step*dt;
					}
				}

				// move apart
				t = (1.0/d) * t;
				double eps = fabs((kulki[i].R + kulki[j].R) - d);
				double dr = eps/2.0;
				r1 = r1 + dr*t;
				r2 = r2 - dr*t;		
				kulki[i].r = r1;
				kulki[j].r = r2;

				// collision
				m1 = kulki[i].R;// 1000.0 * 3.1415 * RAD[i] * RAD[i];
				m2 = kulki[j].R;// 1000.0 * 3.1415 * RAD[j] * RAD[j];
				/*double v1x = kulki[i].v.x;
				double v1y = kulki[i].v.y;
				double v2x = kulki[j].v.x;
				double v2y = kulki[j].v.y;				
				double v1xn = v1x-((2.0 * m2) / (m1 + m2)) * (((v1x - v2x) * (x1 - x2) + (v1y - v2y) * (y1 - y2)) / sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))) * (x1 - x2);
				double v1yn = v1y-((2.0 * m2) / (m1 + m2)) * (((v1x - v2x) * (x1 - x2) + (v1y - v2y) * (y1 - y2)) / sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2))) * (y1 - y2);
				double v2xn = v2x-((2.0 * m1) / (m1 + m2)) * (((v2x - v1x) * (x2 - x1) + (v2y - v1y) * (y2 - y1)) / sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))) * (x2 - x1);
				double v2yn = v2y-((2.0 * m1) / (m1 + m2)) * (((v2x - v1x) * (x2 - x1) + (v2y - v1y) * (y2 - y1)) / sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))) * (y2 - y1);
				kulki[i].v.x = v1xn;
				kulki[i].v.y = v1yn;
				kulki[j].v.x = v2xn;
				kulki[j].v.y = v2yn;		
				*/
				CVector3 v1 = kulki[i].v;
				CVector3 v2 = kulki[j].v;
				kulki[i].v = v1 - (r1-r2)*((2.0 * m2) / (m1 + m2)) * ((v1-v2)*(r1-r2))*(1.0/((r1-r2)*(r1-r2)));
				kulki[j].v = v2 - (r2-r1)*((2.0 * m1) / (m1 + m2)) * ((v2-v1)*(r2-r1))*(1.0/((r2-r1)*(r2-r1)));

				/*CVector3 v1 = kulki[i].v;
				CVector3 v2 = kulki[j].v;
				kulki[i].v = v2;
				kulki[j].v = v1;*/

//				sum2 = kulki[i].v.getLength() + kulki[j].v.getLength();

				/*if(sum1 != sum2)
				{
//					cout << 10000*(sum1 - sum2) << endl;
					double diff = sum2 - sum1;		// vel difference
					CVector3 v1 = kulki[i].v;
					CVector3 v2 = kulki[j].v;
					v1.doNormalize();
					v2.doNormalize();
					kulki[i].v = kulki[i].v + 0.5*diff*v1;
					kulki[j].v = kulki[j].v + 0.5*diff*v2;
				}*/

			}
		}
	}
}


void spread()
{
	for(int i=0; i<kulki.size(); i++)
	if(kulki[i].s!=PASSEDAWAY)
	{
		for(int j=i+1; j<kulki.size(); j++)
		if(kulki[j].s!=PASSEDAWAY)
		{
			float d = (kulki[i].r - kulki[j].r).getLength();
			if(d<R*2)
			if(ran01(mt)<PSPREAD)
			{
				if(kulki[i].s==ILL && kulki[j].s==HEALTHY) 
				{
					kulki[j].s=ILL;
					kulki[j].changet0=step*dt;
				}
				if(kulki[j].s==ILL && kulki[i].s==HEALTHY)
				{
					kulki[i].s=ILL;
					kulki[i].changet0=step*dt;
				}
			}
		}
	}
}
void doconvalescent()
{
	for(int i=0; i<kulki.size(); i++)
	if(kulki[i].s!=PASSEDAWAY)
	{
		if(kulki[i].s==ILL)
		{
			if(step*dt - kulki[i].changet0 > timetoconvalescent)
			{
				if(ran01(mt) < PPASSEDAWAY)
				{
					//kulki[i].s=PASSEDAWAY;
					kulki[i] = kulki[kulki.size()-1];
					kulki.pop_back();
				}
				else
				{
					kulki[i].s=CONVALESCENT;				
					kulki[i].changet0=step*dt;
				}
			}
		}
		if(kulki[i].s==CONVALESCENT)
		{
			if(step*dt - kulki[i].changet0 > convalescenttohealthytime)
				kulki[i].s = HEALTHY;
		}
	}
}
/*		if(kulki[i].s==CONVALESCENT)
		{
			if(step*dt - kulki[i].changet0 > convalescenttohealthytime)
			{
				kulki[i].s=HEALTHY;
			}
		}
		*/


void movelockdown()
{

	for(int k=0; k<N; k++)
	if(ran01(mt)<PLOCKDOWNFLIGHT)
	{


		//int k = int(floor(ran01(mt)*kulki.size()));

		int NPROB = 100;
		int proba = 0;
		int colliding=0;
		float x,y;
		do // find empty spot
		{	
			proba++;
			colliding=0;
			x = ran01(mt);
			y = ran01(mt);
			for(int j=0; j<kulki.size(); j++)
				if((CVector3(x,y,0) - kulki[j].r).getLength() < kulki[j].R )
					colliding = 1;
		} while(colliding==1 && proba < NPROB);
	
		if(proba < NPROB)
		{
			kulki[k].r.x = x;
			kulki[k].r.y = y;
		}	
	}
}


void birth()
{
	if(ran01(mt)>PBIRTH)
		return;

	const float DV = 0.001;
	float vx = ran0505(mt)*DV;
	float vy = ran0505(mt)*DV;
	float _R = R;//+0.1*ran01(mt)*R;	// slight variation
	float x,y;
	int colliding;
	//cout << i << endl;
	int NPROB = 100;
	int proba = 0;
	do // find empty spot
	{	
		proba++;
		colliding=0;
		x = ran01(mt);
		y = ran01(mt);
		for(int j=0; j<kulki.size(); j++)
			if((CVector3(x,y,0) - kulki[j].r).getLength() < _R )
				colliding = 1;
	} while(colliding==1 && proba < NPROB);
	
	if(proba < NPROB)
	{
		kulki.push_back(ball(CVector3(x,y,0),CVector3(vx, vy,0),_R,HEALTHY));
		if(ran01(mt)<PROCSZCZEP)	// zaszczepiony
			kulki[kulki.size()-1].s = VACCINATED;
	}	
}

void count(int &healthy, int &ill, int &convalescent)
{
	ill=0; healthy=0; convalescent=0;
	for(int i=0; i<kulki.size(); i++)
	{
		if(kulki[i].s==HEALTHY) healthy++;
		if(kulki[i].s==VACCINATED) healthy++;
		if(kulki[i].s==ILL) ill++;
		if(kulki[i].s==CONVALESCENT) 
		{
			healthy++;
			convalescent++;
		}
	}
}

void randomill()
{
	int k=kulki.size()*(rand()/float(RAND_MAX));
	kulki[k].s = ILL;
	kulki[k].changet0 = step*dt;
}
void reset()
{
	for(int i=0; i<kulki.size(); i++)
		kulki[i].s=HEALTHY;

	step=0;
}

int findimin(int x, int y)
{
			float dmin = 1e10;
			int imin=0;
			for(int i=0; i<kulki.size(); i++)
			{	
				float d = sqrt(pow((x/(float)W-kulki[i].r.x),2)+pow((y/float(H)-kulki[i].r.y),2));
				if(d<dmin)
				{ 
					dmin=d;
					imin = i;
				}
			}
			if(dmin < R*3)
				return imin;
			else
				return -1;
}



void timerFunction(int data)
{
	step++;
	cout << step << endl;
	if(step==100)
	{
		cout << "NOW" << endl;
		//for(int i=0; i<40; i++)
		//{
			kulki[0].s=ILL;
			kulki[0].changet0=step*dt;
		//}
	}
	// compute
	move();
	collision();
	//birth();
	//spread();
	//doconvalescent();

	int h,i,c;
	static int h0=0,i0=0,c0=0;

	count(h,i,c);
	if(h0!=h || i0!=i || c0!=c)	// only at change
	{
		//file << step << " " << h << " " << i << " " << c << endl;
		cout << step << " " << h << " " << i << " " << c << endl;
	}
	
	glutTimerFunc(0, timerFunction, -1);
}
void key(unsigned char key, int a, int b)
{
	if(int(key)==27) exit(0);
	if(int(key)==114) {init();}		//r
	cout << "key: " << int(key) << endl;
}

void mousemove(int xx, int yy)
{
	
}

void mouse(int button, int state, int xx, int yy)
{
}

int main(int argc, char**argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE| GLUT_RGBA);
	glutInitWindowPosition(30,30);
	cout << "W/H = " << W << "/" << H << endl;
	glutInitWindowSize(W,H);
	glutCreateWindow("Bomb");
	init();
	glutDisplayFunc(render);
	glutKeyboardFunc(key);
    glutMouseFunc(mouse);
	glutMotionFunc( mousemove );
	glutReshapeFunc(redisplay);
	glutTimerFunc(15, timerFunction, -1);
	glutMainLoop();
}

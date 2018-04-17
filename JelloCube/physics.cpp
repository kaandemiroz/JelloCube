/*

USC/Viterbi/Computer Science
"Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"

double restLengthDefault = 0.142857;

#define ABS(x) ((x) < 0 ? -(x) : (x))

/* Computes acceleration to every control point of the jello cube,
which is in state given by 'jello'.
Returns result in array 'a'. */
void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
	int i, j, k, ip, jp, kp;
	double size, dot, restLength, length;
	point p, force, elastic, normal, damping, delta;

#define PROCESS_NEIGHBOUR(di,dj,dk) \
    ip=i+(di);\
    jp=j+(dj);\
    kp=k+(dk);\
    if\
    (!	((ip>7) || (ip<0) || (jp>7) || (jp<0) || (kp>7) || (kp<0)) && \
		((i==0) || (i==7) || (j==0) || (j==7) || (k==0) || (k==7)) && \
		((ip==0) || (ip==7) || (jp==0) || (jp==7) || (kp==0) || (kp==7)) ) \
    {\
		pDIFFERENCE(jello->p[i][j][k], jello->p[ip][jp][kp], normal); \
		pNORMALIZE(normal); \
		pMULTIPLY(normal, length - restLength, elastic); \
		pDIFFERENCE(jello->v[i][j][k], jello->v[ip][jp][kp], damping); \
		pDOT(damping, normal, dot); \
		pMULTIPLY(normal, dot, damping); \
		force.x -= jello->kElastic * elastic.x + jello->dElastic * damping.x; \
		force.y -= jello->kElastic * elastic.y + jello->dElastic * damping.y; \
		force.z -= jello->kElastic * elastic.z + jello->dElastic * damping.z; \
	}\

	// Spring forces
	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				force = { 0, 0, 0 };
				p = jello->p[i][j][k];

				// Structural
				{
					restLength = restLengthDefault;
					PROCESS_NEIGHBOUR(1, 0, 0);
					PROCESS_NEIGHBOUR(0, 1, 0);
					PROCESS_NEIGHBOUR(0, 0, 1);
					PROCESS_NEIGHBOUR(-1, 0, 0);
					PROCESS_NEIGHBOUR(0, -1, 0);
					PROCESS_NEIGHBOUR(0, 0, -1);
				}

				// Shear
				{
					restLength = restLengthDefault * sqrt(2);
					PROCESS_NEIGHBOUR(1, 1, 0);
					PROCESS_NEIGHBOUR(-1, 1, 0);
					PROCESS_NEIGHBOUR(-1, -1, 0);
					PROCESS_NEIGHBOUR(1, -1, 0);
					PROCESS_NEIGHBOUR(0, 1, 1);
					PROCESS_NEIGHBOUR(0, -1, 1);
					PROCESS_NEIGHBOUR(0, -1, -1);
					PROCESS_NEIGHBOUR(0, 1, -1);
					PROCESS_NEIGHBOUR(1, 0, 1);
					PROCESS_NEIGHBOUR(-1, 0, 1);
					PROCESS_NEIGHBOUR(-1, 0, -1);
					PROCESS_NEIGHBOUR(1, 0, -1);

					restLength = restLengthDefault * sqrt(3);
					PROCESS_NEIGHBOUR(1, 1, 1);
					PROCESS_NEIGHBOUR(-1, 1, 1);
					PROCESS_NEIGHBOUR(-1, -1, 1);
					PROCESS_NEIGHBOUR(1, -1, 1);
					PROCESS_NEIGHBOUR(1, 1, -1);
					PROCESS_NEIGHBOUR(-1, 1, -1);
					PROCESS_NEIGHBOUR(-1, -1, -1);
					PROCESS_NEIGHBOUR(1, -1, -1);
				}

				// Bend
				{
					restLength = restLengthDefault * 2;
					PROCESS_NEIGHBOUR(2, 0, 0);
					PROCESS_NEIGHBOUR(0, 2, 0);
					PROCESS_NEIGHBOUR(0, 0, 2);
					PROCESS_NEIGHBOUR(-2, 0, 0);
					PROCESS_NEIGHBOUR(0, -2, 0);
					PROCESS_NEIGHBOUR(0, 0, -2);
				}


				// Bounding Box Collisions
				{
					size = 4;

					if (p.x < -size)
						force.x -= jello->kCollision * p.x - (-size) + jello->dCollision * jello->v[i][j][k].x;
					else if (p.x > size)
						force.x -= jello->kCollision * p.x - (size) + jello->dCollision * jello->v[i][j][k].x;

					if (p.y < -size)
						force.y -= jello->kCollision * p.y - (-size) + jello->dCollision * jello->v[i][j][k].y;
					else if (p.y > size)
						force.y -= jello->kCollision * p.y - (size) + jello->dCollision * jello->v[i][j][k].y;

					if (p.z < -size)
						force.z -= jello->kCollision * p.z - (-size) + jello->dCollision * jello->v[i][j][k].z;
					else if (p.z > size)
						force.z -= jello->kCollision * p.z - (size) + jello->dCollision * jello->v[i][j][k].z;
				}

				// Mouse click force, pushes forward from camera position
				if (g_iLeftMouseButton)
				{
					int power = 40;
					force.x += power * -cos(Phi) * cos(Theta);
					force.y += power * -sin(Phi) * cos(Theta);
					force.z += power * -sin(Theta);
				}

				// Inclined Plane Collisions
				if (jello->incPlanePresent)
				{
					double plane, magnitude;
					normal = { jello->a, jello->b, jello->c };
					pNORMALIZE(normal);
					plane = jello->a * p.x + jello->b * p.y + jello->c * p.z + jello->d;
					if (plane < 0)
					{
						point dist = { p.x + jello->d / jello->a, p.y, p.z};
						pDOT(dist, normal, magnitude);
						pMULTIPLY(normal, magnitude, delta);
						force.x -= jello->kCollision * delta.x + jello->dCollision * jello->v[i][j][k].x;
						force.y -= jello->kCollision * delta.y + jello->dCollision * jello->v[i][j][k].y;
						force.z -= jello->kCollision * delta.z + jello->dCollision * jello->v[i][j][k].z;
					}
				}

				// Force field
				{
					force.x += jello->forceField[i * jello->resolution ^ 2 + j * jello->resolution + k].x;
					force.y += jello->forceField[i * jello->resolution ^ 2 + j * jello->resolution + k].y;
					force.z += jello->forceField[i * jello->resolution ^ 2 + j * jello->resolution + k].z;
				}

				// Final acceleration calculation using F = ma
				a[i][j][k].x = force.x / jello->mass;
				a[i][j][k].y = force.y / jello->mass;
				a[i][j][k].z = force.z / jello->mass;
			}
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
	int i, j, k;
	point a[8][8][8];

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

			}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
	point F1p[8][8][8], F1v[8][8][8],
		F2p[8][8][8], F2v[8][8][8],
		F3p[8][8][8], F3v[8][8][8],
		F4p[8][8][8], F4v[8][8][8];

	point a[8][8][8];


	struct world buffer;

	int i, j, k;

	buffer = *jello; // make a copy of jello

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				pMULTIPLY(jello->v[i][j][k], jello->dt, F1p[i][j][k]);
				pMULTIPLY(a[i][j][k], jello->dt, F1v[i][j][k]);
				pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F2p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F2p[i][j][k]);
				// F2v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F2v[i][j][k]);
				pMULTIPLY(F2p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F2v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F3p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F3v[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);


	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F4p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F4v[i][j][k]);

				pMULTIPLY(F2p[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1p[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4p[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->p[i][j][k], jello->p[i][j][k]);

				pMULTIPLY(F2v[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4v[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->v[i][j][k], jello->v[i][j][k]);
			}

	return;
}

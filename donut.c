#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <float.h>


/* symbols */
#define variousCharacters 12
#define symbolRatio 0.6666666666666666666 /* (8 / 12) - ASCII symbol ratio */
const char illumination[variousCharacters + 1] = ".,-~:;=!*#$@";


/* screen */
#define width 150
#define height 50


/* animation */
#define phiStep 0.006
#define thetaStep 0.06
#define alphaStep 0.006
#define speedMultiplierX 1.
#define speedMultiplierY 2.
#define speedMultiplierZ 3.


/* scene */
#define x_0 (width / 2)    // Shift to draw the torus in the center
#define y_0 (height / 2)   // Shift to draw the torus in the center
#define z_0 35             // Torus Z spawn coordinate. This is literally the point Z=35 in our coordinate system. 
#define scale 8            // Scaling
#define lightVectorX 0
#define lightVectorY 0.7071067811865475 /* (1 / sqrt(2)) */
#define lightVectorZ 0.7071067811865475 /* (1 / sqrt(2)) */
#define torusR 25.0
#define torusr 4.0


/* 
* In the scene, the axes are oriented as follows:
* 
*                            /  Z
*                        |  /  
*                        | /
*                        |/
*           ------------------------------>
*                       /|                X
*                      / |
*                     /  |
*                        V  Y
*                        
*/


void rotateByX(double* x, double* y, double* z, double cos, double sin) {
    double xx = *x, yy = *y, zz = *z;
    *x = xx;
    *y = yy * cos + zz * sin;
    *z = yy * -sin + zz * cos;
} 


void rotateByY(double* x, double* y, double* z, double cos, double sin) {
    double xx = *x, yy = *y, zz = *z;
    *x = xx * cos - zz * sin;
    *y = yy;
    *z = xx * sin + zz * cos;
} 

void rotateByZ(double* x, double* y, double* z, double cos, double sin) {
    double xx = *x, yy = *y, zz = *z;
    *x = xx * cos + yy * sin;
    *y = xx * -sin + yy * cos;
    *z = zz;
} 

void draw_donut(double r1, double r2, double A, double B, double C) {
    /*
    * r1 - Radius of a circle that forms a torus by rotation.
    * r2 - Torus radius. Those this is some X (radius of the hole in this torus) + r1
    *
    * A - rotation angle (radians) along the X axis 
    * B - rotation angle (radians) along the Y axis
    * C - rotation angle (radians) along the Z axis
    */

    const double cosA = cos(A), cosB = cos(B), cosC = cos(C);
    const double sinA = sin(A), sinB = sin(B), sinC = sin(C);
    
    char screen[height][width];
    double zBuff[height][width];

    for (int i = 0; i < height; ++i) {
        memset(screen[i], ' ', width);
        for (int j = 0; j < width; ++j) {
            zBuff[i][j] = DBL_MAX;
        }
    }

    for (double theta = 0; theta < 2 * M_PI; theta += thetaStep) {
        double costheta = cos(theta);
        double sintheta = sin(theta);


        for (double phi = 0; phi < 2 * M_PI; phi += phiStep) {
            double x = r1 * costheta + r2;
            double y = r1 * sintheta;
            double z = 0;
            double nx = costheta;
            double ny = sintheta;
            double nz = 0;

            double sinphi = sin(phi);
            double cosphi = cos(phi);

            /* 
            *  Rotationss to get:
            *     1) A new point on the torus
            *     2) The normal vector to this point 
            */
            rotateByY(&x, &y, &z, cosphi, sinphi);
            rotateByY(&nx, &ny, &nz, cosphi, sinphi);

            /* Rotations for rotation of a point on a torus and its normal */
            rotateByX(&x, &y, &z, cosA, sinA);
            rotateByY(&x, &y, &z, cosB, sinB);
            rotateByZ(&x, &y, &z, cosC, sinC);

            rotateByX(&nx, &ny, &nz, cosA, sinA);
            rotateByY(&nx, &ny, &nz, cosB, sinB);
            rotateByZ(&nx, &ny, &nz, cosC, sinC);

            /* Cosine between the light vector and the normal to the torus */
            double cosalpha = nx * lightVectorX + ny * lightVectorY + nz * lightVectorZ;

            x = x * scale / (z + z_0) / symbolRatio + x_0;
            y = y * scale / (z + z_0) * symbolRatio + y_0;
            
            assert((int)x < width);
            assert((int)y < height);
            assert((int)x >= 0);
            assert((int)y >= 0);
            
            if (cosalpha <= 0 && z < zBuff[((int)y)][(int)x]) {
                zBuff[((int)y)][(int)x] = z;
                screen[((int)y)][(int)x] = illumination[(int)(fabs(cosalpha) * variousCharacters)];
            }
        }
    }

    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            putchar(screen[i][j]);
        }
        putchar('\n');
    }
}

int main() {
    for (double alpha = 0;; alpha += alphaStep) {
        draw_donut(torusr, torusR, alpha * speedMultiplierX, alpha * speedMultiplierY, alpha * speedMultiplierZ);
    }
}

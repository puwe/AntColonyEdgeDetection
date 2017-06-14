/*
 * Copyright 2002-2010 Guillaume Cottenceau.
 *
 * This software may be freely redistributed under the terms
 * of the X11 license.
 *
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#define PNG_DEBUG 3
#include <png.h>
// pheromone
float tau_init = 1e-4;
float alpha = 2;
float beta = 2;
// pheromone evaporation rate
float rho = 0.02;
// construction iterations
unsigned int N = 3;
// theoretical probability
float p_th = 0.01;
// pheromone control thredshold
float t = 0.04;
// block size
unsigned int bx = 256;
unsigned int by = 256;

void abort_(const char * s, ...)
{
        va_list args;
        va_start(args, s);
        vfprintf(stderr, s, args);
        fprintf(stderr, "\n");
        va_end(args);
        abort();
}

unsigned int x, y;

unsigned int Width, Height;
png_byte color_type;
png_byte bit_depth;
png_byte number_of_channels;

png_structp png_ptr;
png_infop info_ptr;
int number_of_passes;
png_bytep * row_pointers;

void read_png_file(char* file_name)
{
        char header[8];    // 8 is the maximum size that can be checked

        /* open file and test for it being a png */
        FILE *fp = fopen(file_name, "rb");
        if (!fp)
                abort_("[read_png_file] File %s could not be opened for reading", file_name);
        fread(header, 1, 8, fp);
        if (png_sig_cmp(header, 0, 8))
                abort_("[read_png_file] File %s is not recognized as a PNG file", file_name);


        /* initialize stuff */
        png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[read_png_file] png_create_read_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[read_png_file] png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during init_io");

        png_init_io(png_ptr, fp);
        png_set_sig_bytes(png_ptr, 8);

        png_read_info(png_ptr, info_ptr);

        Width = png_get_image_width(png_ptr, info_ptr);
        Height = png_get_image_height(png_ptr, info_ptr);
        color_type = png_get_color_type(png_ptr, info_ptr);
        bit_depth = png_get_bit_depth(png_ptr, info_ptr);
        number_of_channels = png_get_channels(png_ptr, info_ptr);

        number_of_passes = png_set_interlace_handling(png_ptr);
        png_read_update_info(png_ptr, info_ptr);


        /* read file */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during read_image");

        row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * Height);
        for (y=0; y<Height; y++)
                row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));

        png_read_image(png_ptr, row_pointers);

        fclose(fp);
}


void write_png_file(char* file_name)
{
        /* create file */
        FILE *fp = fopen(file_name, "wb");
        if (!fp)
                abort_("[write_png_file] File %s could not be opened for writing", file_name);


        /* initialize stuff */
        png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

        if (!png_ptr)
                abort_("[write_png_file] png_create_write_struct failed");

        info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr)
                abort_("[write_png_file] png_create_info_struct failed");

        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during init_io");

        png_init_io(png_ptr, fp);


        /* write header */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing header");

        png_set_IHDR(png_ptr, info_ptr, Width, Height,
                     bit_depth, color_type, PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

        png_write_info(png_ptr, info_ptr);


        /* write bytes */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during writing bytes");

        png_write_image(png_ptr, row_pointers);


        /* end write */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[write_png_file] Error during end of write");

        png_write_end(png_ptr, NULL);

        /* cleanup heap allocation */
        for (y=0; y<Height; y++)
                free(row_pointers[y]);
        free(row_pointers);

        fclose(fp);
}


void ACED_Block(unsigned int x0, unsigned int y0, unsigned int width, unsigned int height, float* Eta)
{
        // image area
        unsigned int A = height*width;
        // image perimeter
        unsigned int P = 2*(height+width);
        if(A==0){
            return;
        }
        
        
        // Initialization
        // ant number
        unsigned int K = floor(sqrt(A));
        float rho_tmp = 1-rho;
        // ant moving steps
        unsigned int L = 3*K;
        // ant memory
        unsigned int M = floor(sqrt(P));
        // K ants moving L steps
        unsigned int ants[K][L];
        //unsigned int* ants = (unsigned int*) malloc(K*L*sizeof(unsigned int));
        //if(NULL==ants){ printf("ants allocation failure!\n"); }
        unsigned int i,j,l,m,k,n,r,q,c;
        int u,v,s;
        unsigned int w;	
        printf("N: %d alpha: %.2f beta: %.2f rho: %.2f t: %.2f p_th: %.2f\n", N, alpha, beta, rho, t, p_th); 
        float tau[height*width];
        //float* tau = (float*) malloc(height*width*sizeof(float));
        //if(NULL==tau){ printf("tau allocation failure!\n"); }
        for(i=0; i<height; i++){
            for(j=0; j<width; j++){
                tau[width*i+j] = tau_init;
            }
        }
        printf("Initialized pheromone matrix tau to %f\n",tau_init);
        float eta[height*width];
        //float* eta = (float*) malloc(height*width*sizeof(float));
        //if(NULL==eta){ printf("eta allocation failure!\n"); }
        unsigned int tc=0;
        for(i=y0; i<y0+height; i++){
            for(j=x0; j<x0+width; j++){
                float e = Eta[i*Width+j];
                if(e>=t){tc++;}
                eta[width*(i-y0)+j-x0] = e;
            }
        }
        if(0==tc){
            if(number_of_channels >= 3){
                for (y=y0; y<y0+height; y++) {
                    png_byte* row = row_pointers[y];
                    for (x=x0; x<x0+width; x++) {
                            png_byte* ptr = &(row[x*number_of_channels]);
                            ptr[0] = 255;
                            ptr[1] = 255;
                            ptr[2] = 255;
                    }
                        
                }
            }else{
                for (y=y0; y<y0+height; y++) {
                    png_byte* row = row_pointers[y];
                    for (x=x0; x<x0+width; x++) {
                            png_byte* ptr = &(row[x*number_of_channels]);
                            ptr[0] = 255;
                    }
                        
                }

            }
            return; 
        }
        for(k=0; k<K; k++){
            c = 0;
            do{
                r = rand()%A;
                for(q=0; q<k; q++){
                    if(ants[q][0] == r || eta[r] < t){
                        break;                       
                    }
                }
                c++;
            }while(q<k&&c<A);
            ants[k][0]=r;
        }
        printf("%d ants with memory length %d and steps %d are randomly dispatched.\n",K,M,L);

        for(n=0; n<N; n++){            
            for(k=0; k<K; k++){
                q = ants[k][0];
                i = q/width;
                j = q%width;
                
                for(l=1; l<L; l++){
                    float p_max=0;
                    float p[4]; // clock-wise
                    if(i>0&&i<height-1&&j>0&&j<width-1){
			u=-1;
			v=0;
			p[0] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=1;
			p[1] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

			u=1;
			v=0;
			p[2] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=-1;
			p[3] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
 
                    }
                    if(i==0&&j>0&&j<width-1){
			u=0;
			v=1;
			p[1] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

			u=1;
			v=0;
			p[2] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=-1;
			p[3] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

                    }
                    if(i==height-1&&j>0&&j<width-1){
			u=-1;
			v=0;
			p[0] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=1;
			p[1] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

			u=0;
			v=-1;
			p[3] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

                    }
                    if(j==0&&i>0&&i<height-1){
			u=-1;
			v=0;
			p[0] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=1;
			p[1] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

			u=1;
			v=0;
			p[2] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
				
                    }
                    if(j==width-1&&i>0&&i<height-1){
			u=-1;
			v=0;
			p[0] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=1;
			v=0;
			p[2] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=-1;
			p[3] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

                    }
                    if(i==0&&j==0){		
			u=0;
			v=1;
			p[1] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

			u=1;
			v=0;
			p[2] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
				
                    }
                    if(i==height-1&&j==0){
			u=-1;
			v=0;
			p[0] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=1;
			p[1] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
                    }
                    if(i==0&&j==width-1){
			u=1;
			v=0;
			p[2] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=-1;
			p[3] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

                    }
                    if(i==height-1&&j==width-1){
			u=-1;
			v=0;
			p[0] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
			
			u=0;
			v=-1;
			p[3] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);

                    }

                    for(m=0; m<M; m++){
                        s = l-1-m;
                        if(s>=0){
                            w = ants[k][s];
                            u = w/width - i;
                            v = w%width - j;
			    if(-1==u&&0==v){
		                p[0] = 0;
			    }
			    if(0==u&&1==v){
				p[1] = 0;
			    }
			    if(1==u&&0==v){
				p[2] = 0;
			    }
			    if(0==u&&-1==v){
				p[3] = 0;
			    } 
                        }else{
                            break;
                        }

                    }
                    unsigned int c_max = 4;
		    for(c=0; c<4; c++){
                        if(p[c]>p_max){
                            p_max = p[c];
                            c_max = c;
                        }
                    }
		    
		    float p_r = (float) rand()/RAND_MAX;
                    if(p_r<p_th){
                        c_max = rand()%4;
                    }

                    if(c_max==4){
                        c = 0;
                        do{
                            r = rand()%A;
                            for(q=0; q<k; q++){
                                if(ants[q][0] == r || eta[r] < t){
                                    break;                       
                                }
                            }
                            c++;
                        }while(q<k&&c<A);
                        ants[k][l] = r;
                        i = r/width;
                        j = r%width;
                        continue;
                    }

                    switch(c_max){
			case 0: u=-1; v=0; break;
			case 1: u=0; v=1; break;
			case 2: u=1; v=0; break;
			case 3: u=0; v=-1; break;
		    }

                    int i_new = i + u;
                    int j_new = j + v;

                    if(i_new>=0&&i_new<height&&j_new>=0&&j_new<width){
                        q = i_new*width + j_new;
                        if(eta[q]>=t){
                            ants[k][l] = q;
                            i = i_new;
                            j = j_new;
                        }
                        else{
                            c = 0;
                            do{
                                r = rand()%A;
                                for(q=0; q<k; q++){
                                    if(ants[q][0] == r || eta[r] < t){
                                        break;                       
                                    }
                                }
                                c++;
                            }while(q<k&&c<A);
                            ants[k][l] = r;
                            i = r/width;
                            j = r%width;
                        }
                    }else{
                        c = 0;
                        do{
                            r = rand()%A;
                            for(q=0; q<k; q++){
                                if(ants[q][0] == r || eta[r] < t){
                                    break;                       
                                }
                            }
                            c++;
                        }while(q<k&&c<A);
                        ants[k][l] = r;
                        i = r/width;
                        j = r%width;

                    }
                }
                        
            }

            for(i=0; i<height; i++){
                for(j=0; j<width; j++){
                    tau[i*width+j] = rho_tmp*tau[i*width+j];
                }
            }

            for(k=0; k<K; k++){
                for(l=0; l<L; l++){
                    q = ants[k][l];
                    tau[q] += eta[q];
                }
                ants[k][0]=q; //ants[k][L-1];
            } 
        }

        // decision
        if(number_of_channels >= 3){
            for (y=y0; y<y0+height; y++) {
                png_byte* row = row_pointers[y];
                for (x=x0; x<x0+width; x++) {
                        png_byte* ptr = &(row[x*number_of_channels]);
                        
                        float g = tau[(y-y0)*width+x-x0];
                        float e = eta[(y-y0)*width+x-x0];
                        if(g>tau_init&&e>=t){
                            ptr[0] = 0;
                            ptr[1] = 0;
                            ptr[2] = 0;
                        }else{
                            ptr[0] = 255;
                            ptr[1] = 255;
                            ptr[2] = 255;
                        }
                                               
                }
                    
            }
        }else{
            for (y=y0; y<y0+height; y++) {
                png_byte* row = row_pointers[y];
                for (x=x0; x<x0+width; x++) {
                        png_byte* ptr = &(row[x*number_of_channels]);
                        
                        float g = tau[(y-y0)*width+x-x0];
                        float e = eta[(y-y0)*width+x-x0];
                        if(g>tau_init&&e>=t){
                            ptr[0] = 0;
                        }else{
                            ptr[0] = 255;
                        }
                                               
                }
                    
            }

        }
        //free(ants);
        //free(tau);
        //free(eta);
}


int main(int argc, char **argv)
{
        if (argc != 3)
                abort_("Usage: program_name <file_in> <file_out>");

        read_png_file(argv[1]);
        
        srand(time(NULL));
        clock_t start = clock(), diff;

        unsigned int y0,x0;
        unsigned int i,j;
        int u,v;
        // Grayscale
        float* IM = (float*) malloc(Height*Width*sizeof(float));
        if(NULL==IM){ printf("IM allocation failure!\n"); }
        float Imax = 0;
        if(number_of_channels >= 3){
            for (y=y0; y<y0+Height; y++) {
                    png_byte* row = row_pointers[y];
                    for (x=x0; x<x0+Width; x++) {
                            png_byte* ptr = &(row[x*number_of_channels]);
                            float g = 0.2989*ptr[0]+0.5870*ptr[1]+0.1140*ptr[2];
                            IM[Width*(y-y0)+(x-x0)] = g;
                            if(g>Imax){
                                Imax = g;
                            }
                            
                    }
            }
        }else{
            for (y=y0; y<y0+Height; y++) {
                    png_byte* row = row_pointers[y];
                    for (x=x0; x<x0+Width; x++) {
                            png_byte* ptr = &(row[x*number_of_channels]);
                            float g = ptr[0];
                            IM[Width*(y-y0)+(x-x0)] = g;
                            if(g>Imax){
                                Imax = g;
                            }
                            
                    }
            }

        }
        printf("Grayscale conversion of %dx%d image finished\n", Height, Width);
        float * Eta = (float*) malloc(Height*Width*sizeof(float));
        if(NULL==Eta){ printf("Eta allocation failure!\n"); }
        float I_diff;
        float I_diff_x, I_diff_y;
	for(i=1; i<Height-1; i++){
            for(j=1; j<Width-1; j++){
                I_diff_x = IM[(i+1)*Width+j]-IM[(i-1)*Width+j];
		I_diff_y = IM[i*Width+j+1]-IM[i*Width+j-1];
                I_diff = 0.5*sqrt(I_diff_x*I_diff_x+I_diff_y*I_diff_y);
                Eta[i*Width+j] = I_diff/Imax;
            }
        }

	i=0;
       	for(j=1; j<Width-1; j++){
	    I_diff_x = IM[(i+1)*Width+j]-IM[i*Width+j];
	    I_diff_y = IM[i*Width+j+1]-IM[i*Width+j-1];
	    I_diff = sqrt(I_diff_x*I_diff_x+0.25*I_diff_y*I_diff_y);
	    Eta[i*Width+j] = I_diff/Imax;
	}
	
	i=Height-1;
	for(j=1; j<Width-1; j++){
	    I_diff_x = IM[i*Width+j]-IM[(i-1)*Width+j];
	    I_diff_y = IM[i*Width+j+1]-IM[i*Width+j-1];
	    I_diff = sqrt(I_diff_x*I_diff_x+0.25*I_diff_y*I_diff_y);
	    Eta[i*Width+j] = I_diff/Imax;
	}
	
	j=0;
	for(i=1; i<Height-1; i++){
	    I_diff_x = IM[(i+1)*Width+j]-IM[(i-1)*Width+j];
	    I_diff_y = IM[i*Width+j+1]-IM[i*Width+j];
	    I_diff = sqrt(0.25*I_diff_x*I_diff_x+I_diff_y*I_diff_y);
	    Eta[i*Width+j] = I_diff/Imax;
        }
	
	j=Width-1;
	for(i=1; i<Height-1; i++){
	    I_diff_x = IM[(i+1)*Width+j]-IM[(i-1)*Width+j];
	    I_diff_y = IM[i*Width+j]-IM[i*Width+j-1];
	    I_diff = sqrt(0.25*I_diff_x*I_diff_x+I_diff_y*I_diff_y);
	    Eta[i*Width+j] = I_diff/Imax;
        }
	
	i=0;
	j=0;
	I_diff_x = IM[(i+1)*Width+j]-IM[i*Width+j];
	I_diff_y = IM[i*Width+j+1]-IM[i*Width+j];
        I_diff = sqrt(I_diff_x*I_diff_x+I_diff_y*I_diff_y);
        Eta[i*Width+j] = I_diff/Imax;
	
	i=0;
	j=Width-1;
	I_diff_x = IM[(i+1)*Width+j]-IM[i*Width+j];
	I_diff_y = IM[i*Width+j]-IM[i*Width+j-1];
        I_diff = sqrt(I_diff_x*I_diff_x+I_diff_y*I_diff_y);
        Eta[i*Width+j] = I_diff/Imax;
	
	i=Height-1;
	j=0;
	I_diff_x = IM[i*Width+j]-IM[(i-1)*Width+j];
	I_diff_y = IM[i*Width+j+1]-IM[i*Width+j];
        I_diff = sqrt(I_diff_x*I_diff_x+I_diff_y*I_diff_y);
        Eta[i*Width+j] = I_diff/Imax;

	i=Height-1;
	j=Width-1;
	I_diff_x = IM[i*Width+j]-IM[(i-1)*Width+j];
	I_diff_y = IM[i*Width+j]-IM[i*Width+j-1];
        I_diff = sqrt(I_diff_x*I_diff_x+I_diff_y*I_diff_y);
        Eta[i*Width+j] = I_diff/Imax;


        printf("Initialized heuristic information eta\n");
        free(IM);

        for(y0=0; y0<Height/by*by; y0+=by){
            for(x0=0; x0<Width/bx*bx; x0+=bx){
                ACED_Block(x0,y0,bx,by,Eta);
                printf("Block 1\n");
            }
            for(x0=Width/bx*bx; x0<Width; x0+=bx){
                ACED_Block(x0,y0,Width-Width/bx*bx,by,Eta);
                printf("Block 2\n");
            }
        }
        for(y0=Height/by*by; y0<Height; y0+=by){
            for(x0=0; x0<Width/bx*bx; x0+=bx){
                ACED_Block(x0,y0,bx,Height-Height/by*by,Eta);
                printf("Block 3\n");
            }
            for(x0=Width/bx*bx; x0<Width; x0+=bx){
                ACED_Block(x0,y0,Width-Width/bx*bx,Height-Height/by*by,Eta);
                printf("Block 4\n");
            }

        }
        //ACED_Block(b,b,b,b);
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
        
        write_png_file(argv[2]);

        return 0;
}

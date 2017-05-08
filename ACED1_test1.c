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

unsigned int width, height;
png_byte color_type;
png_byte bit_depth;

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

        width = png_get_image_width(png_ptr, info_ptr);
        height = png_get_image_height(png_ptr, info_ptr);
        color_type = png_get_color_type(png_ptr, info_ptr);
        bit_depth = png_get_bit_depth(png_ptr, info_ptr);

        number_of_passes = png_set_interlace_handling(png_ptr);
        png_read_update_info(png_ptr, info_ptr);


        /* read file */
        if (setjmp(png_jmpbuf(png_ptr)))
                abort_("[read_png_file] Error during read_image");

        row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
        for (y=0; y<height; y++)
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

        png_set_IHDR(png_ptr, info_ptr, width, height,
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
        for (y=0; y<height; y++)
                free(row_pointers[y]);
        free(row_pointers);

        fclose(fp);
}


void process_file(void)
{
        //if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGB)
        //        abort_("[process_file] input file is PNG_COLOR_TYPE_RGB but must be PNG_COLOR_TYPE_RGBA "
        //               "(lacks the alpha channel)");

        if (png_get_color_type(png_ptr, info_ptr) != PNG_COLOR_TYPE_RGB)
                abort_("[process_file] color_type of input file must be PNG_COLOR_TYPE_RGB (%d) (is %d)",
                       PNG_COLOR_TYPE_RGB, png_get_color_type(png_ptr, info_ptr));

        // Grayscale
        float* IM = (float*) malloc(height*width*sizeof(float));
        float Imax = 0;
        for (y=0; y<height; y++) {
                png_byte* row = row_pointers[y];
                for (x=0; x<width; x++) {
                        png_byte* ptr = &(row[x*3]);
                        float g = 0.2989*ptr[0]+0.5870*ptr[1]+0.1140*ptr[2];
                        //float g = ptr[0];
                        IM[width*y+x] = g;
                        if(g>Imax){
                            Imax = g;
                        }
                }
        }
        printf("Grayscale conversion of %dx%d image finished\n", height, width);
        // Initialization
        // image area
        unsigned long int A = height*width;
        // image perimeter
        unsigned long int P = 2*(height+width);
        // ant number
        unsigned long int K = floor(sqrt(A));
        // pheromone
        float tau_init = 1e-4;
        float alpha = 2;
        float beta = 2;
        // pheromone evaporation rate
        float rho = 0.02;
        // ant moving steps
        unsigned long int L = 3*K;
        // ant memory
        unsigned long int M = floor(sqrt(P));
        // construction iterations
        unsigned int N = 3;
        // theoretical probability
        float p_th = 0.01;
        // pheromone control thredshold
        float t = 0.1;
        // K ants moving L steps
        unsigned long int* ants = (unsigned long int*) malloc(K*L*sizeof(unsigned long int));
        unsigned long int i,j,l,m,k,n,r,q,c;
        printf("N: %d alpha: %.2f beta: %.2f rho: %.2f t: %.2f p_th: %.2f\n", N, alpha, beta, rho, t, p_th); 
        float* tau = (float*) malloc(height*width*sizeof(float));
        for(i=0; i<height; i++){
            for(j=0; j<width; j++){
                tau[width*i+j] = tau_init;
            }
        }
        printf("Initialized pheromone matrix tau to %f\n",tau_init);
        
        float* eta = (float*) malloc(height*width*sizeof(float));
        float I_diff;
        int u,v,w,s;
        for(i=2; i<height-2; i++){
            for(j=2; j<width-2; j++){
                float I_diff_max=0;
                for(u=-2; u<=2; u++){
                    for(v=-2; v<=2; v++){
                        if(u!=0||v!=0){
                            I_diff = fabs(IM[(i+u)*width+(j+v)]-IM[(i-u)*width+(j-v)]);
                            if(I_diff>I_diff_max){
                                I_diff_max = I_diff;
                            }
                        }
                    }
                }
                eta[i*width+j] = I_diff_max/Imax;
            }
        }
        for(i=1; i<height; i += height-3){
            for(j=1; j<width-1; j++){
                float I_diff_max=0;
                for(u=-1; u<=1; u++){
                    for(v=-1; v<=1; v++){
                        if(u!=0||v!=0){
                            I_diff = fabs(IM[(i+u)*width+(j+v)]-IM[(i-u)*width+(j-v)]);
                            if(I_diff>I_diff_max){
                                I_diff_max = I_diff;
                            }
                        }
                    }
                }
                eta[i*width+j] = I_diff_max/Imax;
            }
        }
        for(i=2; i<height-2; i++){
            for(j=1; j<width; j += width-3){
                float I_diff_max=0;
                for(u=-1; u<=1; u++){
                    for(v=-1; v<=1; v++){
                        if(u!=0||v!=0){
                            I_diff = fabs(IM[(i+u)*width+(j+v)]-IM[(i-u)*width+(j-v)]);
                            if(I_diff>I_diff_max){
                                I_diff_max = I_diff;
                            }
                        }
                    }
                }
                eta[i*width+j] = I_diff_max/Imax;
            }
        }
        for(i=0; i<height; i += height-1){
            for(j=1; j<width-1; j++){
                float I_diff_max=0;
                for(v=-1; v<=1&&v!=0; v++){ 
                    I_diff = fabs(IM[i*width+(j+v)]-IM[i*width+(j-v)]);
                    if(I_diff>I_diff_max){
                        I_diff_max = I_diff;
                    }
                }
                eta[i*width+j] = I_diff_max/Imax;
            }
        }
        for(i=1; i<height-1; i++){
            for(j=0; j<width; j += width-1){
                float I_diff_max=0;
                for(u=-1; u<=1&&u!=0; u++){
                    I_diff = fabs(IM[(i+u)*width+j]-IM[(i-u)*width+j]);
                    if(I_diff>I_diff_max){
                        I_diff_max = I_diff;
                    }
                }
                eta[i*width+j] = I_diff_max/Imax;
            }
        } 
        printf("Initialized heuristic information eta\n");

        for(k=0; k<K; k++){
            c = 0;
            do{
                r = rand()%A;
                for(q=0; q<k; q++){
                    if(ants[q*K+0] == r || eta[r] < t){
                        break;                       
                    }
                }
                c++;
            }while(q<k&&c<A);
            ants[k*K+0]=r;
        }
        printf("%d ants with memory length %d and steps %d are randomly dispatched.\n",K,M,L);

        for(n=0; n<N; n++){            
            for(k=0; k<K; k++){
                q = ants[k*K+0];
                i = q/width;
                j = q%width;
                
                for(l=1; l<L; l++){
                    float p_max=0;
                    float p[9];
                    if(i>0&&i<height-1&&j>0&&j<width-1){
                        for(u=-1; u<=1; u++){
                            for(v=-1; v<=1; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }
                    if(i==0&&j>0&&j<width-1){
                        for(u=0; u<=1; u++){
                            for(v=-1; v<=1; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }

                    }
                    if(i==height-1&&j>0&&j<width-1){
                        for(u=-1; u<=0; u++){
                            for(v=-1; v<=1; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }

                    }
                    if(j==0&&i>0&&i<height-1){
                        for(u=-1; u<=1; u++){
                            for(v=0; v<=1; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }
                    if(j==width-1&&i>0&&i<height-1){
                        for(u=-1; u<=1; u++){
                            for(v=-1; v<=0; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }
                    if(i==0&&j==0){
                        for(u=0; u<=1; u++){
                            for(v=0; v<=1; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }
                    if(i==height-1&&j==0){
                        for(u=-1; u<=0; u++){
                            for(v=0; v<=1; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }
                    if(i==0&&j==width-1){
                        for(u=0; u<=1; u++){
                            for(v=-1; v<=0; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }
                    if(i==height-1&&j==width-1){
                        for(u=-1; u<=0; u++){
                            for(v=-1; v<=0; v++){
                                if(u!=0||v!=0){
                                    c = 3*(u+1)+(v+1);
                                    p[c] = pow(tau[(i+u)*width+j+v], alpha)*pow(eta[(i+u)*width+j+v], beta);
                                }
                            }
                        }
                    }

                    for(m=0; m<M; m++){
                        s = l-1-m;
                        if(s>=0){
                            w = ants[k*K+s];
                            u = w/width - i;
                            v = w%width - j;
                            if(1==fmax(abs(u),abs(v))){
                                c = 3*(u+1)+(v+1);
                                p[c] = 0;
                            }
                        }else{
                            break;
                        }

                    }
                    unsigned int c_max = 4;
                    for(c=0; c<9; c++){
                        if(p[c]>p_max){
                            p_max = p[c];
                            c_max = c;
                        }
                    }

                    float p_r = (float) rand()/RAND_MAX;
                    if(p_r<p_th){
                        c_max = rand()%9;
                    }

                    if(c_max==4){
                        c = 0;
                        do{
                            r = rand()%A;
                            for(q=0; q<k; q++){
                                if(ants[q*K+0] == r || eta[r] < t){
                                    break;                       
                                }
                            }
                            c++;
                        }while(q<k&&c<A);
                        ants[k*K+l] = r;
                        i = r/width;
                        j = r%width;
                        continue;
                    }

                    u = c_max/3-1;
                    v = c_max%3-1;

                    long int i_new = i + u;
                    long int j_new = j + v;

                    if(i_new>=0&&i_new<height&&j_new>=0&&j_new<width){
                        q = i_new*width + j_new;
                        if(eta[q]>=t){
                            ants[k*K+l] = q;
                            i = i_new;
                            j = j_new;
                        }
                        else{
                            c = 0;
                            do{
                                r = rand()%A;
                                for(q=0; q<k; q++){
                                    if(ants[q*K+0] == r || eta[r] < t){
                                        break;                       
                                    }
                                }
                                c++;
                            }while(q<k&&c<A);
                            ants[k*K+l] = r;
                            i = r/width;
                            j = r%width;
                        }
                    }else{
                        c = 0;
                        do{
                            r = rand()%A;
                            for(q=0; q<k; q++){
                                if(ants[q*K+0] == r || eta[r] < t){
                                    break;                       
                                }
                            }
                            c++;
                        }while(q<k&&c<A);
                        ants[k*K+l] = r;
                        i = r/width;
                        j = r%width;

                    }
                }
                        
            }

            for(i=0; i<height; i++){
                for(j=0; j<width; j++){
                    tau[i*width+j] = (1-rho)*tau[i*width+j];
                }
            }

            for(k=0; k<K; k++){
                for(l=0; l<L; l++){
                    q = ants[k*K+l];
                    tau[q] += eta[q];
                }
                ants[k*K+0]=ants[k*K+L-1];
            } 
        }

        // decision
        for (y=0; y<height; y++) {
            png_byte* row = row_pointers[y];
            for (x=0; x<width; x++) {
                    png_byte* ptr = &(row[x*3]);
                    
                    float g = tau[y*width+x];
                    float e = eta[y*width+x];
                    if(g>tau_init&&e>t){
                        ptr[0] = 0;
                        ptr[1] = 0;
                        ptr[2] = 0;
                    }else{
                        ptr[0] = 255;
                        ptr[1] = 255;
                        ptr[2] = 255;
                    }
                                           
                    //printf("%f ",g);
            }
                //printf("\n");
                
        }
                
}


int main(int argc, char **argv)
{
        if (argc != 3)
                abort_("Usage: program_name <file_in> <file_out>");

        read_png_file(argv[1]);
        
        srand(time(NULL));

        clock_t start = clock(), diff;
        process_file();
        diff = clock() - start;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("Time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
        
        write_png_file(argv[2]);

        return 0;
}

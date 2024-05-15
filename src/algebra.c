/*
 * @Author: shichenyi shichenyi23@163.com
 * @Date: 2024-05-04 16:23:08
 * @LastEditors: shichenyi shichenyi23@163.com
 * @LastEditTime: 2024-05-15 17:55:51
 * @FilePath: \GitHub\shichenyi.hw1\src\algebra.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "../inc/algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b)
{
    Matrix c;
    int i,j;
    if(a.rows==b.rows&&a.cols==b.cols){
    c = create_matrix(a.rows,a.cols);
    for(i=0;i<a.rows;i++){
        for(j=0;j<a.cols;j++){
            c.data[i][j]=a.data[i][j]+b.data[i][j];
        }
    }
    return c;
    }
    else{
        printf("Error: Matrix a and b must have the same rows and cols.");
        return create_matrix(0, 0);
    }
}

Matrix sub_matrix(Matrix a, Matrix b)
{
    Matrix c;
    int i,j;
    if(a.rows==b.rows&&a.cols==b.cols){
    c = create_matrix(a.rows,a.cols);
    for(i=0;i<a.rows;i++){
        for(j=0;j<a.cols;j++){
            c.data[i][j]=a.data[i][j]-b.data[i][j];
        }
    }
    return c;
    }
    else{
        printf("Error: Matrix a and b must have the same rows and cols.");
        return create_matrix(0, 0);
    }
}

Matrix mul_matrix(Matrix a, Matrix b)
{
    Matrix c;
    int i,j,k;
    if(a.cols==b.rows){
        c = create_matrix(a.rows,b.cols);
        for(i=0;i<a.rows;i++){
            for(j=0;j<b.cols;j++){
                c.data[i][j]=0;
                for(k=0;k<a.cols;k++){
                    c.data[i][j]=c.data[i][j]+a.data[i][k]*b.data[k][j];
                }
            }
        }
        return c;
    }
    else{
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.");
        return create_matrix(0, 0);
    }
}
Matrix scale_matrix(Matrix a, double k)
{
    Matrix c;
    int i,j;
    c = create_matrix(a.rows,a.cols);
    for(i=0;i<a.rows;i++){
        for(j=0;j<a.cols;j++){
            c.data[i][j]=k*a.data[i][j];
        }
    }
    return c;
}

Matrix transpose_matrix(Matrix a)
{
    Matrix c;
    int i,j;
    c = create_matrix(a.cols,a.rows);
    for(i=0;i<a.cols;i++){
        for(j=0;j<a.rows;j++){
            c.data[i][j]=a.data[j][i];
        }
    }
    return c;
}

double det_matrix(Matrix a)
{
   double sum=0;
   int i,j,k;
   Matrix M;
   if(a.rows!=a.cols){
      printf("Error: The matrix must be a square matrix.");
      return 0;
   }
   else{
   if(a.rows==1){
      return a.data[0][0];
   }
   else if(a.rows>1){
      for(i=0;i<a.rows;i++){
        M=create_matrix(a.rows-1,a.cols-1);
        for(j=0;j<a.rows-1;j++){
            for(k=0;k<a.rows-1;k++){
                if(k<i){
                    M.data[j][k]=a.data[j+1][k];
                }
                else if(k>=i){
                    M.data[j][k]=a.data[j+1][k+1];
                }
            }
        }
        sum+=a.data[0][i]*det_matrix(M)*((i%2==1)?-1:1);
      }
      return sum;
   }
   }
}

Matrix inv_matrix(Matrix a)
{
    Matrix Mx,M,Mn;
    int i,j,k,l,m,n;
    double x;
    if(a.rows!=a.cols){
        printf("Error: The matrix must be a square matrix.");
        return create_matrix(0, 0);
    }
    else if(det_matrix(a)==0){
        printf("Error: The matrix is singular.");
        return create_matrix(0, 0);
    }
    else{
        Mx = create_matrix(a.rows,a.cols);
        for(i=0;i<a.rows;i++){
            for(j=0;j<a.rows;j++){
                M = create_matrix(a.rows-1,a.cols-1);
                for(k=0;k<a.rows-1;k++){
                    for(l=0;l<a.rows-1;l++){
                        if(k<i&&l<j){
                            M.data[k][l]=a.data[k][l];
                        }
                        else if(k<i&&l>=j){
                            M.data[k][l]=a.data[k][l+1];
                        }
                        else if(k>=i&&l<j){
                            M.data[k][l]=a.data[k+1][l];
                        }
                        else{
                            M.data[k][l]=a.data[k+1][l+1];
                        }
                    }
                }
                Mx.data[i][j]=det_matrix(M)*(((i+j)%2==1)?-1:1);
            }
        }
        Mn = create_matrix(a.rows,a.cols);
        x = det_matrix(a);
        for(i=0;i<a.rows;i++){
            for(j=0;j<a.rows;j++){
                Mn.data[i][j]=Mx.data[i][j]/x;
            }
        }
        return transpose_matrix(Mn);
    }
}

int rank_matrix(Matrix a)
{
    int i,j,k,l,rank,x;
    double mult,temp;
    rank=a.cols;
    for(i=0;i<rank;i++){
        if(a.data[i][i]!=0){
            for(j=0;j<a.rows;j++){
                if(j!=i){
                    mult=a.data[j][i]/a.data[i][i];
                    for(k=0;k<rank;k++){
                        a.data[j][k]-=mult*a.data[i][k];
                    }
                }
            }
        }
        else{
            x=1;
            for(l=i+1;l<a.rows;l++){
                if(a.data[l][i]!=0){
                    for(j=0;j<rank;j++){
                        temp=a.data[i][j];
                        a.data[i][j]=a.data[l][j];
                        a.data[l][j]=temp;
                    }
                    x=0;
                    break;
                }
            }
            if(x){
                rank--;
                for(l=0;l<a.rows;l++){
                    a.data[l][i]=a.data[l][rank];
                }
            }
            i--;
        }
    }
    return rank;
}
double trace_matrix(Matrix a)
{
    double tr=0;
    int i;
    if(a.rows==a.cols){
        for(i=0;i<a.rows;i++){
            tr=tr+a.data[i][i];
        }
        return tr;
    }
    else{
        printf("Error: The matrix must be a square matrix.");
        return 0;
    }
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}
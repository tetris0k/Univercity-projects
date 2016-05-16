#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include "functions.h"
#include <stdlib.h>
#include <string>
#define EPS 10E-10

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    klava11 = new QShortcut(this);
    klava11->setKey(Qt::Key_1);
    connect(klava11, SIGNAL(activated()), this, SLOT(klava1()));

    klava22 = new QShortcut(this);
    klava22->setKey(Qt::Key_2);
    connect(klava22, SIGNAL(activated()), this, SLOT(klava2()));


    klava33 = new QShortcut(this);
    klava33->setKey(Qt::Key_3);
    connect(klava33, SIGNAL(activated()), this, SLOT(klava3()));

    klava44 = new QShortcut(this);
    klava44->setKey(Qt::Key_4);
    connect(klava44, SIGNAL(activated()), this, SLOT(klava4()));


    klava55 = new QShortcut(this);
    klava55->setKey(Qt::Key_5);
    connect(klava55, SIGNAL(activated()), this, SLOT(klava5()));

    klava66 = new QShortcut(this);
    klava66->setKey(Qt::Key_6);
    connect(klava66, SIGNAL(activated()), this, SLOT(klava6()));

    this->risuy();
}

MainWindow::~MainWindow()
{
    delete ui;
}

bool ff=true;
int met=0;
double grl=-4.;
double grr=4.;
double max=0.;
double min=0.;
int tochg=1000;
int tochin=4;
int tochf=1;
double koef=1.;
std::string file = "input.txt";

void MainWindow::klava1(){
    if (ff) {
        ff=false;
    } else
    {ff=true;}
    MainWindow::risuy();
}

void MainWindow::klava2(){
    met = (met + 1) % 4;
    MainWindow::risuy();
}

void MainWindow::klava3(){
    if (tochin > 1999){
        MainWindow::risuy();
    } else{
    tochin=tochin*2.;
    MainWindow::risuy();
}}

void MainWindow::klava4(){
    if (fabs(koef-1.)<EPS) {koef=0.5;} else if (fabs(koef-0.5)<EPS) {koef=0.25;} else if (fabs(koef-0.25)<EPS) {koef=2.;} else if (fabs(koef-2.)<EPS) {koef=1.5;} else if (fabs(koef-1.5)<EPS) {koef=1.;}
    MainWindow::risuy();
}


void MainWindow::klava5(){
    grl -= 1.;
    MainWindow::risuy();
}

void MainWindow::klava6(){
    grr+=1.;
    MainWindow::risuy();
}

void MainWindow::on_lineEdit_textChanged(const QString &arg1)
{
    bool ok;
    double a = arg1.toDouble(&ok);
    if (ok && grr > a ){
        grl = a;
    }
    MainWindow::risuy();
}

void MainWindow::on_lineEdit_2_textChanged(const QString &arg1)
{
    bool ok;
    double a = arg1.toDouble(&ok);
    if (ok && grl < a ){
        grr = a;
    }
    MainWindow::risuy();
}

void MainWindow::risuy(){
        min=0.;
        max=0.;
        int i=0;
        if (met == 0){
            if(ff){
                QVector<double> x(tochg+2), y(tochg+2);
                double w = (grr - grl)/tochg;
                for (double q= grl; q <= grr; q += w){
                     x[i] = q;
                     y[i] = f(x[i]);
                     if (y[i] > max){
                         max= y[i];
                     }
                     if (y[i] < min){
                         min= y[i];
                     }
                     i++;
                }

                QVector<double> X(tochin), Y(tochin);
                w = (grr-grl)/(tochin-1);
                i = 0;
                for (i = 0; i < tochin; i++){
                     X[i] =grl  + i*w;
                     Y[i] = f(X[i]);
                }
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                QVector<double> a((tochin-1)*100 + 1), b((tochin-1) * 100 + 1);
                MainWindow::metod1(X, Y, a, b, tochin, 1, 0);
                ui->widget->graph(0)->setData(x ,y);
            } else {
                FILE* g;
                g = fopen(file.c_str() , "r");
                if (!g) {
                    ff = true;
                    goto end;
                }
                fscanf(g, "%d", &tochf);
                QVector<double> xf(tochf), yf(tochf);
                for (i = 0; i < tochf; i++){
                    if (!fscanf(g, "%lf", &xf[i])){
                         ff = true;
                         break;
                    }
                    if (!fscanf(g, "%lf", &yf[i])){
                         ff = true;
                         break;
                    }
                    if (yf[i] > max){
                         max= yf[i];
                    }
                    if (yf[i] < min){
                         min= yf[i];
                    }
                }
                fclose (g);
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                QVector<double> a((tochf - 1) * 100 + 1), b((tochf - 1) * 100 + 1);
                MainWindow::metod1(xf, yf, a, b, tochf, 1, 0);
                ui->widget->graph(0)->setData(xf ,yf);

            }
        } else if (met == 1) {
            if (ff){
                QVector<double> x(tochg+2), y(tochg+2);
                double w = (grr - grl)/tochg;
                for (double q = grl; q <= grr; q += w){
                     x[i] = q;
                     y[i] = f(x[i]);
                     if (y[i] > max){
                         max=y[i];
                     }
                     if (y[i] < min){
                         min=y[i];
                     }
                     i++;
                }

                QVector<double> X(tochin), Y(tochin);
                w = (grr - grl)/(tochin-1);
                i = 0;
                for (i = 0; i < tochin; i++){
                     X[i] = grl + i*w;
                     Y[i] = f(X[i]);
                }
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                QVector<double> a((tochin-1)*100+1), b((tochin-1)*100+1);
                MainWindow::metod2(X, Y, a, b, tochin, 1, 0);
                ui->widget->graph(0)->setData(x ,y);
            } else {
                FILE* g;
                g = fopen(file.c_str() , "r");
                if (!g) {
                    ff = true;
                    goto end;
                }
                fscanf(g, "%d", &tochf);
                QVector<double> xf(tochf), yf(tochf);
                for (i = 0; i < tochf; i++){
                    if (!fscanf(g, "%lf", &xf[i])){
                             ff = true;
                             break;
                        }
                        if (!fscanf(g, "%lf", &yf[i])){
                             ff = true;
                             break;
                        }
                        if (yf[i] > max){
                            max = yf[i];
                        }
                        if (yf[i] < min){
                            min= yf[i];
                        }
                    }
                fclose (g);
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                QVector<double> a((tochf - 1) * 100 + 1), b((tochf - 1) * 100 + 1);
                MainWindow::metod2(xf, yf, a, b, tochf, 1, 0);
                ui->widget->graph(0)->setData(xf ,yf);
            }
        } else if (met == 2){
            if(ff){
                QVector<double> x(tochg+2), y(tochg+2);
                double w = (grr-grl)/tochg;
                for (double q = grl;q<=grr;q+=w){
                     x[i] = q;
                     y[i] = f(x[i]);
                     if (y[i]>max){
                         max=y[i];
                     }
                     if (y[i]<min){
                         min=y[i];
                     }
                     i++;
                }
                QVector<double> X(tochin), Y(tochin);
                w = (grr-grl)/(tochin-1);
                i = 0;
                for (i=0;i<tochin;i++){
                     X[i] = grl + i*w;
                     Y[i] = f(X[i]);
                }
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                ui->widget->addGraph();
                ui->widget->graph(0)->setData(x ,y);
                QVector<double> a((tochin-1)*100 + 1), b((tochin-1) * 100 + 1);
                QVector<double> c((tochin-1)*100 + 1), d((tochin-1) * 100 + 1);
                MainWindow::metod1(X, Y, a, b, tochin, 1, 0);
                MainWindow::metod2(X, Y, c, d, tochin, 2, 0);
            } else {
                FILE* g;
                g = fopen(file.c_str() , "r");
                if (!g) {
                    ff = true;
                    goto end;
                }
                fscanf(g, "%d", &tochf);
                QVector<double> xf(tochf), yf(tochf);
                for (i = 0;i < tochf;i++){
                    if (!fscanf(g,"%lf",&xf[i])){
                             ff = true;
                             break;
                        }
                        if (!fscanf(g, "%lf", &yf[i])){
                             ff=true;
                             break;
                        }
                        if (yf[i] > max){
                            max = yf[i];
                        }
                        if (yf[i] < min){
                            min = yf[i];
                        }
                    }
                fclose(g);
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                ui->widget->addGraph();
                ui->widget->graph(0)->setData(xf ,yf);
                QVector<double> a((tochf - 1) * 100 + 1), b((tochf - 1) * 100 + 1);
                QVector<double> c((tochf - 1) * 100 + 1), d((tochf - 1) * 100 + 1);
                MainWindow::metod1(xf, yf, a, b, tochf, 1, 0);
                MainWindow::metod2(xf, yf, c, d, tochf, 2, 0);
            }
        } else {
            if (ff){
                QVector<double> X(tochin), Y(tochin);
                double w= (grr-grl)/(tochin-1);
                i = 0;
                for (i = 0; i<tochin;i++){
                     X[i] = grl + i*w;
                     Y[i] = f(X[i]);
                }
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                QVector<double> a((tochin-1)*100 + 1), b((tochin-1) * 100 + 1);
                QVector<double> c((tochin-1)*100 + 1), d((tochin-1) * 100 + 1);
                MainWindow::metod1(X, Y, a, b, tochin, 0, 1);
                MainWindow::metod2(X, Y, c, d, tochin, 1, 1);
            } else {
                FILE* g;
                g = fopen(file.c_str() , "r");
                if (!g) {
                    ff=true;
                    goto end;
                }
                fscanf(g, "%d", &tochf);
                QVector<double> xf(tochf), yf(tochf);
                for (i = 0; i < tochf; i++){
                    if (!fscanf(g, "%lf", &xf[i])){
                             ff=true;
                             break;
                        }
                        if (!fscanf(g, "%lf", &yf[i])){
                             ff=true;
                             break;
                        }
                    }
                fclose(g);
                ui->widget->clearGraphs();
                ui->widget->addGraph();
                ui->widget->addGraph();
                QVector<double> a((tochf - 1) * 100 + 1), b((tochf - 1) * 100 + 1);
                QVector<double> c((tochf - 1) * 100 + 1), d((tochf - 1) * 100 + 1);
                MainWindow::metod1(xf, yf, a, b, tochf, 0, 1);
                MainWindow::metod2(xf, yf, c, d, tochf, 1, 1);
            }
        }
        ui->widget->xAxis->setLabel("x");
        ui->widget->yAxis->setLabel("y");

        ui->widget->xAxis->setRange(grl * koef, grr * koef);
        ui->widget->yAxis->setRange(koef * min, koef * max);
        ui->widget->replot();
        end:;
        if(ff){
            ui->label_8->setNum(tochin);
        } else {
            ui->label_8->setNum(tochf);
        }
}

double Pi (QVector<double> X, QVector<double> Y, QVector<double> d, int i, double x){
    double a1 = Y[i];
    double a2 = d[i];
    double a3 = ((Y[i+1] - Y[i])/(X[i+1]-X[i]) - a2)/(X[i+1]-X[i]);
    double a4 = (a2 + d[i+1] - 2*(Y[i+1] - Y[i])/(X[i+1]-X[i]))/((X[i+1] - X[i])*(X[i+1] - X[i]));
    return a1 + a2*(x-X[i]) + a3*(x-X[i])*(x-X[i]) + a4*(x-X[i])*(x-X[i])*(x-X[i+1]);
}

double di (QVector<double> x, QVector<double> y, int toch, int i){
    if (i == 0){
        return (3.* fdvoet(x,y,0,1) - di(x,y,toch,1))/2.;
    } else if (i == toch-1){
        return (3.* fdvoet(x,y,toch-2,toch-1) - di(x,y,toch,toch-2))/2.;
    } else {
        return ((x[i+1]-x[i])*fdvoet(x,y,i-1,i) + (x[i]-x[i-1])*fdvoet(x,y,i,i+1))/(x[i+1]-x[i-1]);
    }
}

double fdvoet(QVector<double> x, QVector<double> y, int i, int j){
    return (y[j] - y[i])/(x[j]-x[i]);
}

double Pi_2 (QVector<double> X, QVector<double> Y, double* v, double* E, int i, double x){
    double c1=v[i];
    double c3 = (1./(E[i+1]-X[i]))*((v[i+1]-v[i])/(E[i+1]-E[i]) - (Y[i]-v[i])/(X[i]-E[i]));
    double c2 = (v[i+1]-v[i])/(E[i+1]-E[i]) - c3*(E[i+1]-E[i]);
    return c1+c2*(x-E[i])+c3*(x-E[i])*(x-E[i]);
}

int gde(double* E, int k, double x){
    if (x < E[0]){
        return 0;
    }
    for (int l = 0; l < k-1; l++){
        if ((x-E[l]>-EPS)&&(E[l+1]-x>-EPS)) {
            return l;
        }
    }
    return k-2;
}

void MainWindow::metod1(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int tochf, int kuda, int nev){
    if (ff){
        int i,j,k;
        double w;
        j=0;
        QVector<double> d(tochin);
        for (i=0;i<tochin;i++){
            d[i]=di(X,Y,tochin,i);
        }
        for (i = 0; i < tochin-1; i++){
            w = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*w;
                if (nev == 0){
                    b[j] = Pi(X, Y, d, i, a[j]);
                } else {
                    b[j] = fabs(Pi(X, Y, d, i, a[j])-f(a[j]));
                }
                if (b[j] > max){
                    max= b[j];
                }
                if (b[j] < min){
                    min= b[j];
                }
                j++;
            }
            if (i == tochin-2){
                a[j] = X[tochin-1];
                if (nev==0){
                    b[j] = Pi(X, Y, d, i, a[j]);
                } else {
                    b[j] = fabs(Pi(X, Y, d, i, a[j])-f(a[j]));
                }
                if (b[j] > max){
                    max= b[j];
                }
                if (b[j] <min){
                    min= b[j];
                }
            }
        }
   } else {
        int i,j,k;
        double w;
        j = 0;
        QVector<double> d(tochf);
        for (i = 0; i < tochf; i++){
            d[i] = di(X, Y, tochf, i);
        }
        for (i = 0; i < tochf - 1; i++){
            w = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*w;
                if (nev == 0){
                    b[j] = Pi(X, Y, d, i, a[j]);
                } else {
                    b[j] = fabs(Pi(X, Y, d, i, a[j]) - (Y[i] + k*(Y[i+1]-Y[i])/100.));
                }
                if (b[j] > max){
                    max = b[j];
                }
                if (b[j] < min){
                    min = b[j];
                }
                j++;
            }
            if (i == tochf - 2){
                a[j] = X[tochf - 1];
                if(nev == 0){
                    b[j] = Pi(X, Y, d, i, a[j]);
                } else {
                    b[j] = fabs(Pi(X, Y, d, i, a[j]) - Y[tochf-1]);
                }
                if (b[j] > max){
                    max = b[j];
                }
                if (b[j] < min){
                    min = b[j];
                }
            }
        }
    }
    ui->widget->graph(kuda)->setData(a, b);
    ui->widget->graph(kuda)->setPen(QPen(Qt::red));
}

void MainWindow::metod2(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int tochf, int kuda, int nev){
    if (ff){
        int i,j,k;
        double w;
        j=0;
        double *E = new double [tochin+1];
        for (i=0;i<tochin+1;i++){
            if (i==0){
                E[i]=X[0]-(X[1]-X[0])/2.;
            } else if (i==tochin){
                E[i]=X[tochin-1]+(X[tochin-1]-X[tochin-2])/2.;
            } else {
                E[i]=(X[i-1]+X[i])/2.;
            }
        }
        double *v = new double [tochin+1];
        double *M = new double [(tochin+1)*(tochin+2)];

        double qq=(X[tochin-1]-X[0])/(tochin-1);
        double x0 =X[0]-qq;
        double x00=x0-qq;
        double xl=X[tochin-1]+qq;
        double xll=xl+qq;
        double dx00=(Y[0]-f(x00))/(X[0]-x00);
        double dx01=(Y[2]-Y[0])/(X[2]-X[0]);
        double dxl0=(Y[tochin-1]-Y[tochin-3])/(X[tochin-1]-X[tochin-3]);
        double dxl1=(f(xll)-Y[tochin-1])/(xll-X[tochin-1]);
        double d2fx0=(dx01-dx00)/(X[1]-x0);                      //f''(x[0])
        double d2fxl=(dxl1-dxl0)/(xl-X[tochin-2]);               //f''(x[n-1])

        for (i=0;i<tochin+1;i++){
            if (i==0){
                M[0]=(2./(X[0]-E[0]))/(E[1]-E[0]);
                M[1]=(2./(E[1]-X[0]))/(E[1]-E[0]);
                for (j=2;j<tochin+1;j++){
                    M[j]=0.;
                }
                M[tochin+1]=-1.*(Y[0]*(2./(E[1]-E[0]))*(1./(X[0]-E[0]) + 1./(E[1]-X[0])) + d2fx0);
            } else if (i==tochin){
                for (j=0;j<tochin-1;j++){
                    M[tochin*(tochin+2)+j]=0.;
                }
                M[tochin*(tochin+2)+tochin-1]=(2./(X[tochin-1]-E[tochin-1]))/(E[tochin]-E[tochin-1]);
                M[tochin*(tochin+2)+tochin]=(2./(E[tochin]-X[tochin-1]))/(E[tochin]-E[tochin-1]);
                M[tochin*(tochin+2)+tochin+1]=-1.*(Y[tochin-1]*(2./(E[tochin]-E[tochin-1]))*(1./(X[tochin-1]-E[tochin-1])+1./(E[tochin]-X[tochin-1]))+d2fxl);
            } else {
                for (j = 0; j < i-1; j++){
                    M[i*(tochin+2)+j]=0.;
                }
                M[i*(tochin+2)+i-1]=1./(X[i-1]-E[i-1])-1./(E[i]-E[i-1]);
                M[i*(tochin+2)+i]=1./(E[i]-X[i-1])+1./(E[i]-E[i-1])+1./(X[i]-E[i])+1./(E[i+1]-E[i]);
                M[i*(tochin+2)+i+1]=1./(E[i+1]-X[i])-1./(E[i+1]-E[i]);
                for (j=i+2;j<tochin+1;j++){
                    M[i*(tochin+2)+j]=0.;
                }
                M[i*(tochin+2)+tochin+1]=-1.*(Y[i-1]*(1./(X[i-1]-E[i-1])+1./(E[i]-X[i-1]))+Y[i]*(1./(X[i]-E[i])+1./(E[i+1]-X[i])));
            }
        }
        Solve(tochin+1,M,v);
        for(i=0;i<tochin;i++){
            if(fabs(Pi_2(X,Y,v,E,i,E[i+1]) - v[i+1])>EPS)
                qDebug("very bad\n");
        }
        j=0;
        for (i=0;i<tochin-1;i++){
            if(fabs(Pi_2(X,Y,v,E,i,E[i+1]) - Pi_2(X,Y,v,E,i+1,E[i+1])) > EPS)
                qDebug("bad\n");
            w=(X[i+1]-X[i])/100.;
            for (k=0;k<100;k++){
                a[j] = X[i]+k*w;
                if (nev==0){
                    b[j]=Pi_2(X,Y,v,E,i+(int)(k/50),a[j]);
                } else {
                    b[j]=fabs(Pi_2(X,Y,v,E,gde(E,tochin+1,a[j]),a[j])-f(a[j]));
                }
                if (b[j]>max){
                    max=b[j];
                }
                if (b[j]<min){
                    min=b[j];
                }
                j++;
            }
            if (i==tochin-2){
                a[j]=X[tochin-1];
                if(nev==0){
                    b[j]=Pi_2(X,Y,v,E,tochin-1,a[j]);
                } else {
                    b[j]=fabs(Pi_2(X,Y,v,E,tochin-1,a[j])-f(a[j]));
                }
                if (b[j]>max){
                    max=b[j];
                }
                if (b[j]<min){
                    min=b[j];
                }
            }
        }
        delete []M;
        delete []E;
        delete []v;
    } else {
        int i,j,k;
        double w;
        j = 0;
        double *E = new double [tochf+1];
        for (i = 0; i < tochf+1; i++){
            if (i==0){
                E[i]=X[0]-(X[1]-X[0])/2.;
            } else if (i==tochf){
                E[i]=X[tochf-1]+(X[tochf-1]-X[tochf-2])/2.;
            } else {
                E[i]=(X[i-1]+X[i])/2.;
            }
        }
        double *v = new double [tochf+1];
        double *M = new double [(tochf+1)*(tochf+2)];
        for (i=0;i<tochf+1;i++){
            if (i==0){
                M[0]=(2./(X[0]-E[0]))/(E[1]-E[0]);
                M[1]=(2./(E[1]-X[0]))/(E[1]-E[0]);
                for (j=2;j<tochf+1;j++){
                    M[j]=0.;
                }
                M[tochf+1]=-Y[0]*(2./(E[1]-E[0]))*(1./(X[0]-E[0])+1./(E[1]-X[0]));
            } else if (i==tochf){
                for (j=0;j<tochf-1;j++){
                    M[tochf*(tochf+2)+j]=0.;
                }
                M[tochf*(tochf+2)+tochf-1]=(2./(X[tochf-1]-E[tochf-1]))/(E[tochf]-E[tochf-1]);
                M[tochf*(tochf+2)+tochf]=(2./(E[tochf]-X[tochf-1]))/(E[tochf]-E[tochf-1]);
                M[tochf*(tochf+2)+tochf+1]=-Y[tochf-1]*(2./(E[tochf]-E[tochf-1]))*(1./(X[tochf-1]-E[tochf-1])+1./(E[tochf]-X[tochf-1]));
            } else {
                for (j = 0; j < i-1; j++){
                    M[i*(tochf+2)+j]=0.;
                }
                M[i*(tochf+2)+i-1]=1./(X[i-1]-E[i-1])-1./(E[i]-E[i-1]);
                M[i*(tochf+2)+i]=1./(E[i]-X[i-1])+1./(E[i]-E[i-1])+1./(X[i]-E[i])+1./(E[i+1]-E[i]);
                M[i*(tochf+2)+i+1]=1./(E[i+1]-X[i])-1./(E[i+1]-E[i]);
                for (j=i+2;j<tochf+1;j++){
                    M[i*(tochf+2)+j]=0.;
                }
                M[i*(tochf+2)+tochf+1]=-Y[i-1]*(1./(X[i-1]-E[i-1])+1./(E[i]-X[i-1]))-Y[i]*(1./(X[i]-E[i])+1./(E[i+1]-X[i]));
            }

        }
        Solve(tochf+1,M,v);
        j=0;
        for (i=0;i<tochf-1;i++){
            w=(X[i+1]-X[i])/100.;
            for (k=0;k<100;k++){
                a[j] = X[i]+k*w;
                if (nev==0){
                    b[j]=Pi_2(X,Y,v,E,gde(E,tochin+1,a[j]),a[j]);
                } else {
                    b[j]=fabs(Pi_2(X,Y,v,E,gde(E,tochin+1,a[j]),a[j])-(Y[i]+k*(Y[i+1]-Y[i])/100.));
                }
                if (b[j]>max){
                    max=b[j];
                }
                if (b[j]<min){
                    min=b[j];
                }
                j++;
            }
            if (i==tochf-2){
                a[j]=X[tochf-1];
                if(nev==0){
                    b[j]=Pi_2(X,Y,v,E,gde(E,tochin+1,a[j]),a[j]);
                } else {
                    b[j]=fabs(Pi_2(X,Y,v,E,gde(E,tochin+1,a[j]),a[j])-Y[tochf-1]);
                }
                if (b[j]>max){
                    max=b[j];
                }
                if (b[j]<min){
                    min=b[j];
                }
            }
        }
        delete []M;
        delete []E;
        delete []v;
    }
    ui->widget->graph(kuda)->setData(a, b);
    ui->widget->graph(kuda)->setPen(QPen(Qt::darkGreen));
}

double f(double x){
    return sin(x);
}

void MainWindow::on_lineEdit_3_textChanged(const QString &arg1)
{
    if (arg1 == ""){
        file = "input.txt";
    } else {
        file = arg1.toStdString();
    }
}

void MainWindow::on_pushButton_2_clicked()
{
    ff = false;
    MainWindow::risuy();
}

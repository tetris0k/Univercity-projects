#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    key1 = new QShortcut(this);
    key1->setKey(Qt::Key_1);
    connect(key1, SIGNAL(activated()), this, SLOT(slotShortcut1()));

    key2 = new QShortcut(this);
    key2->setKey(Qt::Key_2);
    connect(key2, SIGNAL(activated()), this, SLOT(slotShortcut2()));


    key3 = new QShortcut(this);
    key3->setKey(Qt::Key_3);
    connect(key3, SIGNAL(activated()), this, SLOT(slotShortcut3()));

    key4 = new QShortcut(this);
    key4->setKey(Qt::Key_4);
    connect(key4, SIGNAL(activated()), this, SLOT(slotShortcut4()));

    this->DrawPlot();
}

MainWindow::~MainWindow()
{
    delete ui;
}

static int mode = 0;
static int func = 0;
static double coef = 1.;
static double left_b = -5.;
static double right_b = 5.;
static int N = 6;
static int n = 1000;
static int n_f = 1;
static double min_func = 0.;
static double max_func = 0.;

void MainWindow::on_lineEdit_textChanged(const QString &arg1)
{
    bool ok;
    double a = arg1.toDouble(&ok);
    if (ok && right_b > a ){
        left_b = a;
    }
    MainWindow::DrawPlot();
}

void MainWindow::on_lineEdit_2_textChanged(const QString &arg1)
{
    bool ok;
    double a = arg1.toDouble(&ok);
    if (ok && left_b < a){
        right_b = a;
    }
    MainWindow::DrawPlot();
}

void MainWindow::on_lineEdit_3_textChanged(const QString &arg1)
{
    bool ok;
    int a = arg1.toInt(&ok, 10);
    if ((ok) && ((mode == 0 && a > 1) || (mode > 0 && a > 3))){
        N = a;
    }
    MainWindow::DrawPlot();
}

void MainWindow::DrawPlot()
{
    max_func = 0.;
    min_func = 0.;
    int i = 0;
    if (func == 0){
        if (mode == 0){
            QVector<double> x(n+2), y(n+2);
            double h = (right_b - left_b)/n;
            for (double tmp = left_b; tmp <= right_b; tmp += h){
                 x[i] = tmp;
                 y[i] = function(x[i]);
                 if (y[i] > max_func){
                     max_func = y[i];
                 }
                 if (y[i] < min_func){
                     min_func = y[i];
                 }
                 i++;
            }

            QVector<double> X(N), Y(N);
            h = (right_b - left_b)/(N-1);
            i = 0;
            for (i = 0; i < N; i++){
                 X[i] = left_b + i*h;
                 Y[i] = function(X[i]);
            }
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            QVector<double> a((N-1)*100 + 1), b((N-1) * 100 + 1);
            MainWindow::FirstMethod(X, Y, a, b, N, 1);
            ui->widget->graph(0)->setData(x ,y);
            ui->widget->graph(0)->setName("Initial function");
            ui->widget->graph(0)->setPen(QPen(Qt::blue));
        } else if (mode == 1){
            QVector<double> x(n+2), y(n+2);
            double h = (right_b - left_b)/n;
            for (double tmp = left_b; tmp <= right_b; tmp += h){
                 x[i] = tmp;
                 y[i] = function(x[i]);
                 if (y[i] > max_func){
                     max_func = y[i];
                 }
                 if (y[i] < min_func){
                     min_func = y[i];
                 }
                 i++;
            }

            QVector<double> X(N), Y(N);
            h = (right_b - left_b)/(N-1);
            i = 0;
            for (i = 0; i < N; i++){
                 X[i] = left_b + i*h;
                 Y[i] = function(X[i]);
            }
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            QVector<double> a((N-1)*100 + 1), b((N-1) * 100 + 1);
            MainWindow::SecondMethod(X, Y, a, b, N, 1);
            ui->widget->graph(0)->setData(x ,y);
            ui->widget->graph(0)->setName("Initial function");
            ui->widget->graph(0)->setPen(QPen(Qt::blue));
        } else if (mode == 2){
            QVector<double> x(n+2), y(n+2);
            double h = (right_b - left_b)/n;
            for (double tmp = left_b; tmp <= right_b; tmp += h){
                 x[i] = tmp;
                 y[i] = function(x[i]);
                 if (y[i] > max_func){
                     max_func = y[i];
                 }
                 if (y[i] < min_func){
                     min_func = y[i];
                 }
                 i++;
            }
            QVector<double> X(N), Y(N);
            h = (right_b - left_b)/(N-1);
            i = 0;
            for (i = 0; i < N; i++){
                 X[i] = left_b + i*h;
                 Y[i] = function(X[i]);
            }
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            ui->widget->addGraph();
            ui->widget->graph(0)->setData(x ,y);
            ui->widget->graph(0)->setName("Initial function");
            ui->widget->graph(0)->setPen(QPen(Qt::blue));
            QVector<double> a((N-1)*100 + 1), b((N-1) * 100 + 1);
            QVector<double> a1((N-1)*100 + 1), b1((N-1) * 100 + 1);
            MainWindow::FirstMethod(X, Y, a, b, N, 1);
            MainWindow::SecondMethod(X, Y, a1, b1, N, 2);
        } else {
            QVector<double> X(N), Y(N);
            double h = (right_b - left_b)/(N-1);
            i = 0;
            for (i = 0; i < N; i++){
                 X[i] = left_b + i*h;
                 Y[i] = function(X[i]);
            }
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            QVector<double> a((N-1)*100 + 1), b((N-1) * 100 + 1);
            QVector<double> a1((N-1)*100 + 1), b1((N-1) * 100 + 1);
            MainWindow::FirstResid(X, Y, a, b, N, 0);
            MainWindow::SecondResid(X, Y, a1, b1, N, 1);
        }
        ui->widget->xAxis->setLabel("x");
        ui->widget->yAxis->setLabel("y");

        ui->widget->xAxis->setRange(left_b * coef, right_b * coef);
        ui->widget->yAxis->setRange(coef * (min_func - 1.), coef * (max_func + 1.));
        ui->widget->replot();
    } else {
        FILE* f;
        f = fopen("input.txt" , "r");
        if (!f) {
            func = 0;
            goto end;
        }
        fscanf(f, "%d", &n_f);
        if (mode == 0){
            QVector<double> x_f(n_f), y_f(n_f);
            for (i = 0; i < n_f; i++){
                if (!fscanf(f, "%lf", &x_f[i])){
                         func = 0;
                         break;
                    }
                    if (!fscanf(f, "%lf", &y_f[i])){
                         func = 0;
                         break;
                    }
                    if (y_f[i] > max_func){
                        max_func = y_f[i];
                    }
                    if (y_f[i] < min_func){
                        min_func = y_f[i];
                    }
                }
            fclose (f);
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            QVector<double> a((n_f - 1) * 100 + 1), b((n_f - 1) * 100 + 1);
            MainWindow::FirstMethod(x_f, y_f, a, b, n_f, 1);
            ui->widget->graph(0)->setData(x_f ,y_f);
            ui->widget->graph(0)->setName("Initial function");
            ui->widget->graph(0)->setPen(QPen(Qt::blue));
        } else if (mode == 1){
            QVector<double> x_f(n_f), y_f(n_f);
            for (i = 0; i < n_f; i++){
                if (!fscanf(f, "%lf", &x_f[i])){
                         func = 0;
                         break;
                    }
                    if (!fscanf(f, "%lf", &y_f[i])){
                         func = 0;
                         break;
                    }
                    if (y_f[i] > max_func){
                        max_func = y_f[i];
                    }
                    if (y_f[i] < min_func){
                        min_func = y_f[i];
                    }
                }
            fclose (f);
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            QVector<double> a((n_f - 1) * 100 + 1), b((n_f - 1) * 100 + 1);
            MainWindow::SecondMethod(x_f, y_f, a, b, n_f, 1);
            ui->widget->graph(0)->setData(x_f ,y_f);
            ui->widget->graph(0)->setName("Initial function");
            ui->widget->graph(0)->setPen(QPen(Qt::blue));
        } else if (mode == 2){
            QVector<double> x_f(n_f), y_f(n_f);
            for (i = 0; i < n_f; i++){
                if (!fscanf(f, "%lf", &x_f[i])){
                         func = 0;
                         break;
                    }
                    if (!fscanf(f, "%lf", &y_f[i])){
                         func = 0;
                         break;
                    }
                    if (y_f[i] > max_func){
                        max_func = y_f[i];
                    }
                    if (y_f[i] < min_func){
                        min_func = y_f[i];
                    }
                }
            fclose (f);
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            ui->widget->addGraph();
            ui->widget->graph(0)->setData(x_f ,y_f);
            ui->widget->graph(0)->setName("Initial function");
            ui->widget->graph(0)->setPen(QPen(Qt::blue));
            QVector<double> a((n_f - 1) * 100 + 1), b((n_f - 1) * 100 + 1);
            QVector<double> a1((n_f - 1) * 100 + 1), b1((n_f - 1) * 100 + 1);
            MainWindow::FirstMethod(x_f, y_f, a, b, n_f, 1);
            MainWindow::SecondMethod(x_f, y_f, a1, b1, n_f, 2);
        } else {
            QVector<double> x_f(n_f), y_f(n_f);
            for (i = 0; i < n_f; i++){
                if (!fscanf(f, "%lf", &x_f[i])){
                         func = 0;
                         break;
                    }
                    if (!fscanf(f, "%lf", &y_f[i])){
                         func = 0;
                         break;
                    }
                    if (y_f[i] > max_func){
                        max_func = y_f[i];
                    }
                    if (y_f[i] < min_func){
                        min_func = y_f[i];
                    }
                }
            fclose (f);
            ui->widget->clearGraphs();
            ui->widget->addGraph();
            ui->widget->addGraph();
            QVector<double> a((n_f - 1) * 100 + 1), b((n_f - 1) * 100 + 1);
            QVector<double> a1((n_f - 1) * 100 + 1), b1((n_f - 1) * 100 + 1);
            MainWindow::FirstResid(x_f, y_f, a, b, n_f, 0);
            MainWindow::SecondResid(x_f, y_f, a1, b1, n_f, 1);
        }
        ui->widget->xAxis->setLabel("x");
        ui->widget->yAxis->setLabel("y");

        ui->widget->xAxis->setRange(left_b * coef, right_b * coef);
        ui->widget->yAxis->setRange(coef * min_func, coef * max_func);
        ui->widget->replot();
    }
    end:;
}

void MainWindow::on_comboBox_activated(int index)
{
    func = index;
    MainWindow::DrawPlot();
}

void MainWindow::on_comboBox_2_activated(int index)
{
    mode = index;
    if (N < 4 && mode > 0){
        N = 4;
        ui->lineEdit_3->setText("4");
    }
    MainWindow::DrawPlot();
}

void MainWindow::on_comboBox_3_activated(int index)
{
    if (index == 0){
        coef = 1.;
    } else if (index == 1){
        coef = 0.75;
    } else if (index == 2){
        coef = 0.5;
    } else if (index == 3){
        coef = 1.5;
    } else {
        coef = 2.;
    }
    MainWindow::DrawPlot();
}

void MainWindow::FirstMethod(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest){
    if (func == 0){
        int i,j,k;
        double h;
        double x0 = X[0] - ((X[N-1] - X[0])/(N-1));
        double fx0 = function(x0);
        double xl = X[N-1] + ((X[N-1] - X[0])/(N-1));
        double fl = function(xl);
        j = 0;
        QVector<double> d(N);
        for (i = 0; i < N; i++){
            d[i] = di(X,Y,i,x0,fx0,xl,fl);
        }
        for (i = 0; i < N-1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = Pi(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == N-2){
                a[j] = X[N-1];
                b[j] = Pi(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
    } else {
        int i,j,k;
        double h;
        j = 0;
        QVector<double> d(n_f);
        for (i = 0; i < n_f; i++){
            d[i] = di_file(X, Y, n_f, i);
        }
        for (i = 0; i < n_f - 1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = Pi(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == n_f - 2){
                a[j] = X[n_f - 1];
                b[j] = Pi(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
    }
    ui->widget->graph(dest)->setData(a, b);
    ui->widget->graph(dest)->setName("First method");
    ui->widget->graph(dest)->setPen(QPen(Qt::red));
}

void MainWindow::SecondMethod(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest){
    if (func == 0){
        int i,j,k;
        double h;
        j = 0;
        double *M = new double [N*(N+1)];
        for (i = 0; i < N; i++){
            if (i == 0){
                M[0] = X[2] - X[1];
                M[1] = X[2] - X[0];
                for (j = 2; j < N; j++){
                    M[j] = 0.;
                }
                M[N] = -1.*((Y[1]-Y[0])*(X[2]-X[1])*(2.*X[2]+X[1]-3.*X[0])/(X[1]-X[0]) + (Y[2]-Y[1])*(X[1]-X[0])*(X[1]-X[0])/(X[2]-X[1]))/(X[2]-X[0]);
            } else if (i == N-1){
                for (j = 0; j < N-2; j++){
                    M[(N-1)*(N+1) + j] = 0.;
                }
                M[N*(N+1) - 3] = X[N-1] - X[N-3];
                M[N*(N+1) - 2] = X[N-2] - X[N-3];
                M[N*(N+1) - 1] = -1.*((Y[N-2] - Y[N-3])*(X[N-1]-X[N-2])*(X[N-1]-X[N-2])/(X[N-2]-X[N-3]) + (Y[N-1]-Y[N-2])*(X[N-2]-X[N-3])*(3.*X[N-1]-X[N-2]-2.*X[N-3])/(X[N-1]-X[N-2]))/(X[N-1] - X[N-3]);
            } else {
                for (j = 0; j < i-1; j++){
                    M[i*(N+1) + j] = 0.;
                }
                M[i*(N+1) + i-1] = X[i+1] - X[i];
                M[i*(N+1) + i] = 2.*(X[i+1]-X[i-1]);
                M[i*(N+1) + i+1] = X[i] - X[i-1];
                for (j = i+2; j < N; j++){
                    M[i*(N+1) + j] = 0.;
                }
                M[i*(N+1) + N] = -3.*(Y[i]-Y[i-1])*(X[i+1]-X[i])/(X[i]-X[i-1]) - 3.*(Y[i+1]-Y[i])*(X[i]-X[i-1])/(X[i+1]-X[i]);
            }
        }
        double* d = new double [N];
        Solve(N,M,d);
        delete []M;
        j = 0;
        for (i = 0; i < N-1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = Pi_2(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == N-2){
                a[j] = X[N-1];
                b[j] = Pi_2(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
        delete []d;
    } else {
        int i,j,k;
        double h;
        j = 0;
        double *M = new double [n_f*(n_f+1)];
        for (i = 0; i < n_f; i++){
            if (i == 0){
                M[0] = X[2] - X[1];
                M[1] = X[2] - X[0];
                for (j = 2; j < n_f; j++){
                    M[j] = 0.;
                }
                M[n_f] = -1.*((Y[1]-Y[0])*(X[2]-X[1])*(2.*X[2]+X[1]-3.*X[0])/(X[1]-X[0]) + (Y[2]-Y[1])*(X[1]-X[0])*(X[1]-X[0])/(X[2]-X[1]))/(X[2]-X[0]);
            } else if (i == n_f-1){
                for (j = 0; j < n_f-2; j++){
                    M[(n_f-1)*(n_f+1) + j] = 0.;
                }
                M[n_f*(n_f+1) - 3] = X[n_f-1] - X[n_f-3];
                M[n_f*(n_f+1) - 2] = X[n_f-2] - X[n_f-3];
                M[n_f*(n_f+1) - 1] = -1.*((Y[n_f-2] - Y[n_f-3])*(X[n_f-1]-X[n_f-2])*(X[n_f-1]-X[n_f-2])/(X[n_f-2]-X[n_f-3]) + (Y[n_f-1]-Y[n_f-2])*(X[n_f-2]-X[n_f-3])*(3.*X[n_f-1]-X[n_f-2]-2.*X[n_f-3])/(X[n_f-1]-X[n_f-2]))/(X[n_f-1] - X[n_f-3]);
            } else {
                for (j = 0; j < i-1; j++){
                    M[i*(n_f+1) + j] = 0.;
                }
                M[i*(n_f+1) + i-1] = X[i+1] - X[i];
                M[i*(n_f+1) + i] = 2.*(X[i+1]-X[i-1]);
                M[i*(n_f+1) + i+1] = X[i] - X[i-1];
                for (j = i+2; j < n_f; j++){
                    M[i*(n_f+1) + j] = 0.;
                }
                M[i*(n_f+1) + n_f] = -3.*(Y[i]-Y[i-1])*(X[i+1]-X[i])/(X[i]-X[i-1]) - 3.*(Y[i+1]-Y[i])*(X[i]-X[i-1])/(X[i+1]-X[i]);
            }
        }
        double* d = new double [n_f];
        Solve(n_f, M, d);
        delete []M;
        j = 0;
        for (i = 0; i < n_f-1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = Pi_2(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == n_f-2){
                a[j] = X[n_f-1];
                b[j] = Pi_2(X, Y, d, i, a[j]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
        delete []d;
    }
    ui->widget->graph(dest)->setData(a, b);
    ui->widget->graph(dest)->setName("Second method");
    ui->widget->graph(dest)->setPen(QPen(Qt::darkGreen));
}

void MainWindow::FirstResid(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest){
    if (func == 0){
        int i,j,k;
        double h;
        double x0 = X[0] - ((X[N-1] - X[0])/(N-1));
        double fx0 = function(x0);
        double xl = X[N-1] + ((X[N-1] - X[0])/(N-1));
        double fl = function(xl);
        j = 0;
        QVector<double> d(N);
        for (i = 0; i < N; i++){
            d[i] = di(X,Y,i,x0,fx0,xl,fl);
        }
        for (i = 0; i < N-1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = fabs(Pi(X, Y, d, i, a[j]) - function(a[j]));
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == N-2){
                a[j] = X[N-1];
                b[j] = fabs(Pi(X, Y, d, i, a[j]) - Y[N-1]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
    } else {
        int i,j,k;
        double h;
        j = 0;
        QVector<double> d(n_f);
        for (i = 0; i < n_f; i++){
            d[i] = di_file(X, Y, n_f, i);
        }
        for (i = 0; i < n_f - 1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = fabs(Pi(X, Y, d, i, a[j]) - (Y[i] + k*(Y[i+1]-Y[i])/100.));
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == n_f - 2){
                a[j] = X[n_f - 1];
                b[j] = fabs(Pi(X, Y, d, i, a[j]) - Y[n_f - 1]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
    }
    ui->widget->graph(dest)->setData(a, b);
    ui->widget->graph(dest)->setName("First residual norm");
    ui->widget->graph(dest)->setPen(QPen(Qt::red));
}

void MainWindow::SecondResid(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest){
    if (func == 0){
        int i,j,k;
        double h;
        j = 0;
        double *M = new double [N*(N+1)];
        for (i = 0; i < N; i++){
            if (i == 0){
                M[0] = X[2] - X[1];
                M[1] = X[2] - X[0];
                for (j = 2; j < N; j++){
                    M[j] = 0.;
                }
                M[N] = -1.*((Y[1]-Y[0])*(X[2]-X[1])*(2.*X[2]+X[1]-3.*X[0])/(X[1]-X[0]) + (Y[2]-Y[1])*(X[1]-X[0])*(X[1]-X[0])/(X[2]-X[1]))/(X[2]-X[0]);
            } else if (i == N-1){
                for (j = 0; j < N-2; j++){
                    M[(N-1)*(N+1) + j] = 0.;
                }
                M[N*(N+1) - 3] = X[N-1] - X[N-3];
                M[N*(N+1) - 2] = X[N-2] - X[N-3];
                M[N*(N+1) - 1] = -1.*((Y[N-2] - Y[N-3])*(X[N-1]-X[N-2])*(X[N-1]-X[N-2])/(X[N-2]-X[N-3]) + (Y[N-1]-Y[N-2])*(X[N-2]-X[N-3])*(3.*X[N-1]-X[N-2]-2.*X[N-3])/(X[N-1]-X[N-2]))/(X[N-1] - X[N-3]);
            } else {
                for (j = 0; j < i-1; j++){
                    M[i*(N+1) + j] = 0.;
                }
                M[i*(N+1) + i-1] = X[i+1] - X[i];
                M[i*(N+1) + i] = 2.*(X[i+1]-X[i-1]);
                M[i*(N+1) + i+1] = X[i] - X[i-1];
                for (j = i+2; j < N; j++){
                    M[i*(N+1) + j] = 0.;
                }
                M[i*(N+1) + N] = -3.*(Y[i]-Y[i-1])*(X[i+1]-X[i])/(X[i]-X[i-1]) - 3.*(Y[i+1]-Y[i])*(X[i]-X[i-1])/(X[i+1]-X[i]);
            }
        }
        double* d = new double [N];
        Solve(N,M,d);
        delete []M;
        j = 0;
        for (i = 0; i < N-1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = fabs(Pi_2(X, Y, d, i, a[j]) - function(a[j]));
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == N-2){
                a[j] = X[N-1];
                b[j] = fabs(Pi_2(X, Y, d, i, a[j]) - Y[N-1]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
        delete []d;
    } else {
        int i,j,k;
        double h;
        j = 0;
        double *M = new double [n_f*(n_f+1)];
        for (i = 0; i < n_f; i++){
            if (i == 0){
                M[0] = X[2] - X[1];
                M[1] = X[2] - X[0];
                for (j = 2; j < n_f; j++){
                    M[j] = 0.;
                }
                M[n_f] = -1.*((Y[1]-Y[0])*(X[2]-X[1])*(2.*X[2]+X[1]-3.*X[0])/(X[1]-X[0]) + (Y[2]-Y[1])*(X[1]-X[0])*(X[1]-X[0])/(X[2]-X[1]))/(X[2]-X[0]);
            } else if (i == n_f-1){
                for (j = 0; j < n_f-2; j++){
                    M[(n_f-1)*(n_f+1) + j] = 0.;
                }
                M[n_f*(n_f+1) - 3] = X[n_f-1] - X[n_f-3];
                M[n_f*(n_f+1) - 2] = X[n_f-2] - X[n_f-3];
                M[n_f*(n_f+1) - 1] = -1.*((Y[n_f-2] - Y[n_f-3])*(X[n_f-1]-X[n_f-2])*(X[n_f-1]-X[n_f-2])/(X[n_f-2]-X[n_f-3]) + (Y[n_f-1]-Y[n_f-2])*(X[n_f-2]-X[n_f-3])*(3.*X[n_f-1]-X[n_f-2]-2.*X[n_f-3])/(X[n_f-1]-X[n_f-2]))/(X[n_f-1] - X[n_f-3]);
            } else {
                for (j = 0; j < i-1; j++){
                    M[i*(n_f+1) + j] = 0.;
                }
                M[i*(n_f+1) + i-1] = X[i+1] - X[i];
                M[i*(n_f+1) + i] = 2.*(X[i+1]-X[i-1]);
                M[i*(n_f+1) + i+1] = X[i] - X[i-1];
                for (j = i+2; j < n_f; j++){
                    M[i*(n_f+1) + j] = 0.;
                }
                M[i*(n_f+1) + n_f] = -3.*(Y[i]-Y[i-1])*(X[i+1]-X[i])/(X[i]-X[i-1]) - 3.*(Y[i+1]-Y[i])*(X[i]-X[i-1])/(X[i+1]-X[i]);
            }
        }
        double* d = new double [n_f];
        Solve(n_f, M, d);
        delete []M;
        j = 0;
        for (i = 0; i < n_f-1; i++){
            h = (X[i+1] - X[i])/100.;
            for (k = 0; k < 100; k++){
                a[j] = X[i] + k*h;
                b[j] = fabs(Pi_2(X, Y, d, i, a[j]) - (Y[i] + k*(Y[i+1]-Y[i])/100.));
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
                j++;
            }
            if (i == n_f-2){
                a[j] = X[n_f-1];
                b[j] = fabs(Pi_2(X, Y, d, i, a[j]) - Y[n_f-1]);
                if (b[j] > max_func){
                    max_func = b[j];
                }
                if (b[j] < min_func){
                    min_func = b[j];
                }
            }
        }
        delete []d;
    }
    ui->widget->graph(dest)->setData(a, b);
    ui->widget->graph(dest)->setName("Second residual norm");
    ui->widget->graph(dest)->setPen(QPen(Qt::darkGreen));
}

double function(double x){
    return x*sin(x);
}

double di (QVector<double> x, QVector<double> y, int i, double x0, double fx0, double xl, double fl){
    if (i == 0){
        return (y[0] - fx0)/(2.*(x[0] - x0)) + (y[1] - y[0])/(2.*(x[1] - x[0]));
    } else if (i == N-1){
        return (y[N-1] - y[N-2])/(2.*(x[N-1]-x[N-2])) + (fl - y[N-1])/(2.*(xl - x[N-1]));
    } else {
        return (y[i] - y[i-1])/(2.*(x[i]-x[i-1])) + (y[i+1] - y[i])/(2.*(x[i+1] - x[i]));
    }
}

double di_file (QVector<double> x, QVector<double> y, int n_f, int i){
    if (i == 0){
        return (y[1] - y[0])/(x[1] - x[0]);
    } else if (i == n_f - 1){
        return (y[n_f - 1] - y[n_f - 2])/(x[n_f - 1] - x[n_f - 2]);
    } else {
        return (y[i] - y[i-1])/(2.*(x[i]-x[i-1])) + (y[i+1] - y[i])/(2.*(x[i+1] - x[i]));
    }
}

double Pi (QVector<double> X, QVector<double> Y, QVector<double> d, int i, double x){
    double a1 = Y[i];
    double a2 = d[i];
    double a3 = ((Y[i+1] - Y[i])/(X[i+1]-X[i]) - a2)/(X[i+1]-X[i]);
    double a4 = (a2 + d[i+1] - 2*(Y[i+1] - Y[i])/(X[i+1]-X[i]))/((X[i+1] - X[i])*(X[i+1] - X[i]));
    return a1 + a2*(x-X[i]) + a3*(x-X[i])*(x-X[i]) + a4*(x-X[i])*(x-X[i])*(x-X[i+1]);
}

double Pi_2 (QVector<double> X, QVector<double> Y, double* d, int i, double x){
    double a1 = Y[i];
    double a2 = d[i];
    double a3 = ((Y[i+1] - Y[i])/(X[i+1]-X[i]) - a2)/(X[i+1]-X[i]);
    double a4 = (a2 + d[i+1] - 2*(Y[i+1] - Y[i])/(X[i+1]-X[i]))/((X[i+1] - X[i])*(X[i+1] - X[i]));
    return a1 + a2*(x-X[i]) + a3*(x-X[i])*(x-X[i]) + a4*(x-X[i])*(x-X[i])*(x-X[i+1]);
}

void MainWindow::slotShortcut1(){
    func = (func + 1) % 2;
    ui->comboBox->setCurrentIndex(func);
    MainWindow::DrawPlot();
}

void MainWindow::slotShortcut2(){
    mode = (mode + 1) % 4;
    ui->comboBox_2->setCurrentIndex(mode);
    MainWindow::DrawPlot();
}

void MainWindow::slotShortcut3(){
    int tmp = ui->comboBox_3->currentIndex();
    ui->comboBox_3->setCurrentIndex((tmp + 1) % 5);
    MainWindow::on_comboBox_3_activated((tmp + 1) % 5);
}

void MainWindow::slotShortcut4(){
    if ((mode == 0 && N > 0.5) || (mode > 0 && N > 1)){
        N *= 2;
        ui->lineEdit_3->setText(QString::number(N));
    }
    MainWindow::DrawPlot();
}

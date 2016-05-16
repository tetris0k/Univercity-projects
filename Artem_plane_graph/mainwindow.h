#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QShortcut>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:

    void on_lineEdit_textChanged(const QString &arg1);

    void on_lineEdit_2_textChanged(const QString &arg1);

    void on_lineEdit_3_textChanged(const QString &arg1);

    void on_comboBox_activated(int index);

    void on_comboBox_2_activated(int index);

    void DrawPlot();

    void on_comboBox_3_activated(int index);

    void slotShortcut1();

    void slotShortcut2();

    void slotShortcut3();

    void slotShortcut4();

    void FirstMethod(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest);

    void SecondMethod(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest);

    void FirstResid(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest);

    void SecondResid(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int n_f, int dest);

private:
    Ui::MainWindow *ui;
    QShortcut *key1;
    QShortcut *key2;
    QShortcut *key3;
    QShortcut *key4;
};

double function(double x);
double di (QVector<double> x, QVector<double> y, int i, double x0, double fx0, double xl, double fl);
double di_file (QVector<double> x, QVector<double> y, int n_f, int i);
double Pi (QVector<double> X, QVector<double> Y, QVector<double> d, int i, double x);
double Pi_2 (QVector<double> X, QVector<double> Y, double* d, int i, double x);
double smult(int ,int ,int ,double* , double*);
double smult2(int ,double * , int , double *);
void master(int , double *, int  , int , double *, double *);
void nowv(int ,double *, double *, int );
void Solve(int , double* , double* );
#endif // MAINWINDOW_H

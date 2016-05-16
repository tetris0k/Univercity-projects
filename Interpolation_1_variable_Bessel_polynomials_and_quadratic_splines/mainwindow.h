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

    void risuy();

    void klava1();

    void klava2();

    void klava3();

    void klava4();

    void klava5();

    void klava6();

    void metod1(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int tochf, int kuda, int nev);

    void metod2(QVector<double> X, QVector<double> Y, QVector<double> a, QVector<double> b, int tochf, int kuda, int nev);

    void on_lineEdit_textChanged(const QString &arg1);

    void on_lineEdit_2_textChanged(const QString &arg1);

    void on_lineEdit_3_textChanged(const QString &arg1);

    void on_pushButton_2_clicked();

private:
    Ui::MainWindow *ui;

    QShortcut *klava11;
    QShortcut *klava22;
    QShortcut *klava33;
    QShortcut *klava44;
    QShortcut *klava55;
    QShortcut *klava66;
};

double f(double x);
double di (QVector<double> x, QVector<double> y, int toch, int i);
double Pi (QVector<double> X, QVector<double> Y, QVector<double> d, int i, double x);
double fdvoet(QVector<double> x, QVector<double> y, int i, int j);
double Pi_2 (QVector<double> X, QVector<double> Y, double* v, double* E, int i, double x);
int gde(double* E, int k, double x);
#endif // MAINWINDOW_H


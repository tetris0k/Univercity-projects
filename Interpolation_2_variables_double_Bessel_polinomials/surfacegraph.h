#ifndef SURFACEGRAPH_H
#define SURFACEGRAPH_H

#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtWidgets/QLabel>
#include <QtWidgets/QRadioButton>
#include <QShortcut>

using namespace QtDataVisualization;

class SurfaceGraph : public QObject
{
    Q_OBJECT
public:
    explicit SurfaceGraph(Q3DSurface *surface);
    ~SurfaceGraph();

    QShortcut *key1;
    QShortcut *key2;
    QShortcut *key3;
    QShortcut *key4;
    QShortcut *key5;
    QShortcut *key6;
    QShortcut *key7;
    QShortcut *key8;
    QShortcut *key9;
    int nn;
    int mm;
    int kol;

    void enableInitPlotModel(bool enable);
    void enableInterpModel(bool enable);
    void enableInitInterpModel(bool enable);
    void enableResidModel(bool enable);

    void adjustXMin(int min);
    void adjustXMax(int max);
    void adjustZMin(int min);
    void adjustZMax(int max);

    void setCount(QLabel *countN, QLabel* countM);
    void setKol(QRadioButton *first, QRadioButton *second,QRadioButton *third,QRadioButton *fourth);  //добавляем функцию выставления точек

private slots:
    void slotShortcut1();

    void slotShortcut2();

    void slotShortcut3();

    void slotShortcut4();

    void slotShortcut5();

    void slotShortcut6();

    void slotShortcut7();

    void slotShortcut8();

    void slotShortcut9();

private:
    Q3DSurface *m_graph;
    QSurfaceDataProxy *m_InterpProxy;
    QSurfaceDataProxy *m_InitPlotProxy;
    QSurface3DSeries *m_InterpSeries;
    QSurface3DSeries *m_InitPlotSeries;
    QSurfaceDataProxy *m_ResidProxy;
    QSurface3DSeries *m_ResidSeries;

    QLabel *m_countN;
    QLabel *m_countM;
    double m_rangeMinX;
    double m_rangeMinZ;
    double m_stepX;
    double m_stepZ;

    QRadioButton *m_first;         //добавила несколько кнопок внутрь класса SurfaceGraph
    QRadioButton *m_second;
    QRadioButton *m_third;
    QRadioButton *m_fourth;


    void setAxisXRange(double min, double max);
    void setAxisZRange(double min, double max);
    void fillInitPlotProxy();
    void fillInterpProxy();
    void fillResidProxy();
};

double func(double , double );
double proizvX(double *x, double* f, int i, int j);
double proizvZ(double *z, double* f, int i, int j);
void Aminus(double* A, double h);
void transp(double* a, int u);
void UmnozhMatric(double* a, double* b, double* c);
void Fij(double* f, double* dx, double* dz, double* dxdz, double* F, int i, int j);
double Pf(double**** g, double* x, double* z, int i, int j, double X, double Z);
int gde(double X, double* x, int l);
#endif // SURFACEGRAPH_H

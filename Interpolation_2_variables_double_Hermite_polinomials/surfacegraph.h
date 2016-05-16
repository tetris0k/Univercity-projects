#ifndef SURFACEGRAPH_H
#define SURFACEGRAPH_H

#include <QtDataVisualization/Q3DSurface>
#include <QtDataVisualization/QSurfaceDataProxy>
#include <QtDataVisualization/QSurface3DSeries>
#include <QtWidgets/QSlider>
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

    int nn;
    int mm;

    void enableInitPlotModel(bool enable);
    void enableInterpModel(bool enable);
    void enableInitInterpModel(bool enable);
    void enableResidModel(bool enable);

    void adjustXMin(int min);
    void adjustXMax(int max);
    void adjustZMin(int min);
    void adjustZMax(int max);

    void setAxisMinSliderX(QSlider *slider) { m_axisMinSliderX = slider; }
    void setAxisMaxSliderX(QSlider *slider) { m_axisMaxSliderX = slider; }
    void setAxisMinSliderZ(QSlider *slider) { m_axisMinSliderZ = slider; }
    void setAxisMaxSliderZ(QSlider *slider) { m_axisMaxSliderZ = slider; }

    void setCount(QLabel *countN, QLabel* countM);
    void setMode(QRadioButton *init, QRadioButton *interp, QRadioButton *together, QRadioButton *resid);

private slots:
    void slotShortcut1();

    void slotShortcut2();

    void slotShortcut3();

    void slotShortcut4();

    void slotShortcut5();

private:
    Q3DSurface *m_graph;
    QSurfaceDataProxy *m_InterpProxy;
    QSurfaceDataProxy *m_InitPlotProxy;
    QSurface3DSeries *m_InterpSeries;
    QSurface3DSeries *m_InitPlotSeries;
    QSurfaceDataProxy *m_ResidProxy;
    QSurface3DSeries *m_ResidSeries;

    QRadioButton *m_init;
    QRadioButton *m_interp;
    QRadioButton *m_together;
    QRadioButton *m_resid;
    QLabel *m_countN;
    QLabel *m_countM;
    QSlider *m_axisMinSliderX;
    QSlider *m_axisMaxSliderX;
    QSlider *m_axisMinSliderZ;
    QSlider *m_axisMaxSliderZ;
    double m_rangeMinX;
    double m_rangeMinZ;
    double m_stepX;
    double m_stepZ;


    void setAxisXRange(double min, double max);
    void setAxisZRange(double min, double max);
    void fillInitPlotProxy();
    void fillInterpProxy();
    void fillResidProxy();
};

double function(double , double );
double didx (double *x, double* f, int i, int j, double x0, double fx0, double xl, double fl);
double didz (double *z, double* f, int i, int j, double z0, double fz0, double zl, double fzl);
void Aminus(double* A, double h);
void transp (double* a, int u);
void MatrMult(double* a, double* b, double* c);
void Fij(double* f, double* dx, double* dz, double* dxdz, double* F, int i, int j);
double Pf(double**** g, double* x, double* z, int i, int j, double X, double Z);
int which_in(double X, double* x, int l);
//double d2f(double x, double z, double x_step, double z_step);
#endif // SURFACEGRAPH_H

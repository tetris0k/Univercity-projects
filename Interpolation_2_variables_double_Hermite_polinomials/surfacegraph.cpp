#include "surfacegraph.h"

#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/Q3DTheme>
#include <QtGui/QImage>
#include <QtCore/qmath.h>

#define EPS 10e-10

using namespace QtDataVisualization;

static int CountX = 50;
static int CountZ = 50;
double Min = -2.0;
double Max = 2.0;
int n = 4;
int m = 4;
int how_much = 0;

SurfaceGraph::SurfaceGraph(Q3DSurface *surface)
    : m_graph(surface)
{
    m_graph->setAxisX(new QValue3DAxis);
    m_graph->setAxisY(new QValue3DAxis);
    m_graph->setAxisZ(new QValue3DAxis);

    this->nn = n;
    this->mm = m;

    m_InitPlotProxy = new QSurfaceDataProxy();
    m_InitPlotSeries = new QSurface3DSeries(m_InitPlotProxy);
    fillInitPlotProxy();

    m_InterpProxy = new QSurfaceDataProxy();
    m_InterpSeries = new QSurface3DSeries(m_InterpProxy);
    fillInterpProxy();

    m_ResidProxy = new QSurfaceDataProxy();
    m_ResidSeries = new QSurface3DSeries(m_ResidProxy);
    fillResidProxy();

}

SurfaceGraph::~SurfaceGraph()
{
    delete m_graph;
}

void SurfaceGraph::fillInitPlotProxy()
{
    double stepX = (Max - Min) / double(CountX - 1);
    double stepZ = (Max - Min) / double(CountZ - 1);

    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(CountZ);
    for (int i = 0 ; i < CountZ ; i++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(CountX);
        double z = qMin(Max, (i * stepZ + Min));
        int index = 0;
        for (int j = 0; j < CountX; j++) {
            double x = qMin(Max, (j * stepX + Min));
            double y = function(x,z);
            (*newRow)[index++].setPosition(QVector3D(x, y, z));
        }
        *dataArray << newRow;
    }

    m_InitPlotProxy->resetArray(dataArray);
}

void SurfaceGraph::fillInterpProxy()
{
    double X,Y,Z;
    int i,j,index, k, l;
    double *f = new double [n*m];
    double *x = new double [n];
    double *z = new double [m];
    double intstepx = (Max - Min) / double(n - 1);
    double intstepz = (Max - Min) / double(m - 1);
    for (i = 0; i < n; i++){
        x[i] = Min + i*intstepx;
        for (j = 0; j < m; j++){
            z[j] = Min + j * intstepz;
            f[i*m +j] = function(x[i], z[j]);
        }
    }

    double *ddz = new double [n*m];
    double z0 = z[0] - intstepz;
    double zl = z[m-1] + intstepz;
    double fz0, fzl;
    for (i = 0; i < n; i++){
        fz0 = function(x[i], z0);
        fzl = function(x[i], zl);
        for (j = 0; j < m; j++){
            ddz[i*m + j] = didz(z, f, i, j, z0, fz0, zl, fzl);
        }
    }
    double x0, fx0, xl, fxl;
    x0 = x[0] - intstepx;
    xl = x[n-1] + intstepx;
    double *ddx = new double [n*m];
    for (j = 0; j < m; j++){
        fx0 = function(x0, z[j]);
        fxl = function(xl, z[j]);
        for (i = 0; i < n; i++){
            ddx[i*m + j] = didx(x,f,i,j,x0,fx0,xl,fxl);
        }
    }

    double *ddxdz = new double [n*m];
    double dz0, dzl;
    for (i = 0; i < n; i++){
        if (i == 0){
            dz0 = (function(x[0], z0) - function(x0, z0))/(2.*(x[0] - x0)) + (function(x[1],z0) - function(x[0], z0))/(2.*(x[1] - x[0]));
            dzl = (function(x[0], zl) - function(x0, zl))/(2.*(x[0] - x0)) + (function(x[1],zl) - function(x[0], zl))/(2.*(x[1] - x[0]));
        } else if (i == n-1){
            dz0 = (function(x[n-1], z0) - function(x[n-2], z0))/(2.*(x[n-1] - x[n-2])) + (function(xl,z0) - function(x[n-1], z0))/(2.*(xl - x[n-1]));
            dzl = (function(x[n-1], zl) - function(x[n-2], zl))/(2.*(x[n-1] - x[n-2])) + (function(xl,zl) - function(x[n-1], zl))/(2.*(xl - x[n-1]));
        } else {
            dz0 = (function(x[i], z0) - function(x[i-1], z0))/(2.*(x[i] - x[i-1])) + (function(x[i+1],z0) - function(x[i], z0))/(2.*(x[i+1] - x[i]));
            dzl = (function(x[i], zl) - function(x[i-1], zl))/(2.*(x[i] - x[i-1])) + (function(x[i+1],zl) - function(x[i], zl))/(2.*(x[i+1] - x[i]));
        }
        for (j = 0; j < m; j++){
            ddxdz[i*m + j] = didz(z, ddx, i, j, z0, dz0, zl, dzl);
        }
    }

    double ****g = new double *** [n-1];
    for(i = 0; i < n-1; i++){
        g[i] = new double ** [m-1];
        for (j = 0; j < m-1; j++){
            g[i][j] = new double * [4];
            for (k = 0; k < 4; k++){
                g[i][j][k] = new double [4];
            }
        }
    }
    double *A = new double [16];
    double *B = new double [16];
    double *C = new double [16];
    double *F = new double [16];

    for (i = 0; i < n-1; i++){
        for (j = 0; j < m-1; j++){
            Fij(f,ddx,ddz,ddxdz,F,i,j);
            Aminus(A, fabs(x[i+1] - x[i]));
            Aminus(B, fabs(z[j+1] - z[j]));
            transp(B, 4);
            MatrMult(A, F, C);
            MatrMult(C, B, F);
            for (k = 0; k < 4; k++){
                for (l = 0; l < 4; l++){
                    g[i][j][k][l] = F[k*4 + l];
                }
            }
        }
    }
    delete [] f;
    delete [] ddx;
    delete [] ddz;
    delete [] ddxdz;
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] F;
    int ii, jj;
    double stepX = (Max - Min) / double(CountX - 1);
    double stepZ = (Max - Min) / double(CountZ - 1);
    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(CountZ);
    for (j = 0 ; j < CountZ ; j++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(CountX);
        Z = qMin(Max, (j * stepZ + Min));
        jj = which_in(Z,z,m);
        index = 0;
        for (i = 0; i < CountX; i++) {
            X = qMin(Max, (i * stepX + Min));
            ii = which_in(X,x,n);
            Y = Pf(g, x, z, ii, jj, X, Z);
            (*newRow)[index++].setPosition(QVector3D(X, Y, Z));
        }
        *dataArray << newRow;
    }
    delete [] x;
    delete [] z;
    for(i = 0; i < n-1; i++){
        for (j = 0; j < m-1; j++){
            for (k = 0; k < 4; k++){
                delete [] g[i][j][k];
            }
            delete [] g[i][j];
        }
        delete [] g[i];
    }
    delete [] g;
    m_InterpProxy->resetArray(dataArray);
}

void SurfaceGraph::fillResidProxy(){
    this->nn = n;
    this->mm = m;
    double X,Y,Z;
    int i,j,index, k, l;
    double *f = new double [n*m];
    double *x = new double [n];
    double *z = new double [m];
    double intstepx = (Max - Min) / double(n - 1);
    double intstepz = (Max - Min) / double(m - 1);
    for (i = 0; i < n; i++){
        x[i] = Min + i*intstepx;
        for (j = 0; j < m; j++){
            z[j] = Min + j * intstepz;
            f[i*m +j] = function(x[i], z[j]);
        }
    }
    double x0, fx0, xl, fxl;
    x0 = x[0] - intstepx;
    xl = x[n-1] + intstepx;
    double *ddx = new double [n*m];
    for (j = 0; j < m; j++){
        fx0 = function(x0, z[j]);
        fxl = function(xl, z[j]);
        for (i = 0; i < n; i++){
            ddx[i*m + j] = didx(x,f,i,j,x0,fx0,xl,fxl);
        }
    }
    double *ddz = new double [n*m];
    double z0 = z[0] - intstepz;
    double zl = z[m-1] + intstepz;
    double fz0, fzl;
    for (i = 0; i < n; i++){
        fz0 = function(x[i] , z0);
        fzl = function(x[i], zl);
        for (j = 0; j < m; j++){
            ddz[i*m + j] = didz(z, f, i, j, z0, fz0, zl, fzl);
        }
    }
    double *ddxdz = new double [n*m];
    double dz0, dzl;
    for (i = 0; i < n; i++){
        if (i == 0){
            dz0 = (function(x[0], z0) - function(x0, z0))/(2.*(x[0] - x0)) + (function(x[1],z0) - function(x[0], z0))/(2.*(x[1] - x[0]));
            dzl = (function(x[0], zl) - function(x0, zl))/(2.*(x[0] - x0)) + (function(x[1],zl) - function(x[0], zl))/(2.*(x[1] - x[0]));
        } else if (i == n-1){
            dz0 = (function(x[n-1], z0) - function(x[n-2], z0))/(2.*(x[n-1] - x[n-2])) + (function(xl,z0) - function(x[n-1], z0))/(2.*(xl - x[n-1]));
            dzl = (function(x[n-1], zl) - function(x[n-2], zl))/(2.*(x[n-1] - x[n-2])) + (function(xl,zl) - function(x[n-1], zl))/(2.*(xl - x[n-1]));
        } else {
            dz0 = (function(x[i], z0) - function(x[i-1], z0))/(2.*(x[i] - x[i-1])) + (function(x[i+1],z0) - function(x[i], z0))/(2.*(x[i+1] - x[i]));
            dzl = (function(x[i], zl) - function(x[i-1], zl))/(2.*(x[i] - x[i-1])) + (function(x[i+1],zl) - function(x[i], zl))/(2.*(x[i+1] - x[i]));
        }
        for (j = 0; j < m; j++){
            ddxdz[i*m + j] = didz(z, ddx, i, j, z0, dz0, zl, dzl);
        }
    }

    double ****g = new double *** [n-1];
    for(i = 0; i < n-1; i++){
        g[i] = new double ** [m-1];
        for (j = 0; j < m-1; j++){
            g[i][j] = new double * [4];
            for (k = 0; k < 4; k++){
                g[i][j][k] = new double [4];
            }
        }
    }
    double *A = new double [16];
    double *B = new double [16];
    double *C = new double [16];
    double *F = new double [16];

    for (i = 0; i < n-1; i++){
        for (j = 0; j < m-1; j++){
            Fij(f,ddx,ddz,ddxdz,F,i,j);
            Aminus(A, fabs(x[i+1] - x[i]));
            Aminus(B, fabs(z[j+1] - z[j]));
            transp(B, 4);
            MatrMult(A, F, C);
            MatrMult(C, B, F);
            for (k = 0; k < 4; k++){
                for (l = 0; l < 4; l++){
                    g[i][j][k][l] = F[k*4 + l];
                }
            }
        }
    }
    delete [] f;
    delete [] ddx;
    delete [] ddz;
    delete [] ddxdz;
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] F;
    int ii, jj;
    double stepX = (Max - Min) / double(CountX - 1);
    double stepZ = (Max - Min) / double(CountZ - 1);
    QSurfaceDataArray *dataArray = new QSurfaceDataArray;
    dataArray->reserve(CountZ);
    for (j = 0 ; j < CountZ ; j++) {
        QSurfaceDataRow *newRow = new QSurfaceDataRow(CountX);
        Z = qMin(Max, (j * stepZ + Min));
        jj = which_in(Z,z,m);
        index = 0;
        for (i = 0; i < CountX; i++) {
            X = qMin(Max, (i * stepX + Min));
            ii = which_in(X,x,n);
            Y = fabs(Pf(g, x, z, ii, jj, X, Z) - function(X,Z));
            (*newRow)[index++].setPosition(QVector3D(X, Y, Z));
        }
        *dataArray << newRow;
    }
    delete [] x;
    delete [] z;
    for(i = 0; i < n-1; i++){
        for (j = 0; j < m-1; j++){
            for (k = 0; k < 4; k++){
                delete [] g[i][j][k];
            }
            delete [] g[i][j];
        }
        delete [] g[i];
    }
    delete [] g;
    m_ResidProxy->resetArray(dataArray);
}

int which_in(double X, double* x, int l){
    if (X < x[0]){
        return 0;
    }
    for (int k = 0; k < l-1; k++){
        if (( X - x[k] > -EPS) && ( x[k+1] - X > -EPS )) {
            return k;
        }
    }
    return l-2;
}

void SurfaceGraph::enableInitPlotModel(bool enable)
{
    if (enable) {
        this->nn = n;
        this->mm = m;
        m_InitPlotSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InitPlotSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(Min, Max);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(Min, Max);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);

        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_InitPlotSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemeArmyBlue);

        m_rangeMinX = Min;
        m_rangeMinZ = Min;
        m_stepX = (Max - Min) / double(CountX - 1);
        m_stepZ = (Max - Min) / double(CountZ - 1);
        m_axisMinSliderX->setMaximum(CountX - 2);
        m_axisMinSliderX->setValue(0);
        m_axisMaxSliderX->setMaximum(CountX - 1);
        m_axisMaxSliderX->setValue(CountX - 1);
        m_axisMinSliderZ->setMaximum(CountZ - 2);
        m_axisMinSliderZ->setValue(0);
        m_axisMaxSliderZ->setMaximum(CountZ - 1);
        m_axisMaxSliderZ->setValue(CountZ - 1);
        how_much = 0;
    }
}

void SurfaceGraph::enableInterpModel(bool enable)
{
    if (enable) {
        this->nn = n;
        this->mm = m;
        m_InterpSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InterpSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(Min, Max);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(Min, Max);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);

        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_InterpSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemeArmyBlue);


        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::red);
        m_graph->seriesList().at(0)->setBaseGradient(gr);
        m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

        m_rangeMinX = Min;
        m_rangeMinZ = Min;
        m_stepX = (Max - Min) / double(CountX - 1);
        m_stepZ = (Max - Min) / double(CountZ - 1);
        m_axisMinSliderX->setMaximum(CountX - 2);
        m_axisMinSliderX->setValue(0);
        m_axisMaxSliderX->setMaximum(CountX - 1);
        m_axisMaxSliderX->setValue(CountX - 1);
        m_axisMinSliderZ->setMaximum(CountZ - 2);
        m_axisMinSliderZ->setValue(0);
        m_axisMaxSliderZ->setMaximum(CountZ - 1);
        m_axisMaxSliderZ->setValue(CountZ - 1);
        how_much = 1;
    }
}

void SurfaceGraph::enableInitInterpModel(bool enable)
{
    if (enable) {
        this->nn = n;
        this->mm = m;
        m_InterpSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InterpSeries->setFlatShadingEnabled(true);
        m_InitPlotSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InitPlotSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(Min, Max);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(Min, Max);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);
        m_graph->setSelectionMode(QAbstract3DGraph::SelectionMultiSeries);

        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_InitPlotSeries);
        m_graph->addSeries(m_InterpSeries);

        m_graph->activeTheme()->setType(Q3DTheme::ThemeArmyBlue);

        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::red);
        m_graph->seriesList().at(1)->setBaseGradient(gr);
        m_graph->seriesList().at(1)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

        m_rangeMinX = Min;
        m_rangeMinZ = Min;
        m_stepX = (Max - Min) / double(CountX - 1);
        m_stepZ = (Max - Min) / double(CountZ - 1);
        m_axisMinSliderX->setMaximum(CountX - 2);
        m_axisMinSliderX->setValue(0);
        m_axisMaxSliderX->setMaximum(CountX - 1);
        m_axisMaxSliderX->setValue(CountX - 1);
        m_axisMinSliderZ->setMaximum(CountZ - 2);
        m_axisMinSliderZ->setValue(0);
        m_axisMaxSliderZ->setMaximum(CountZ - 1);
        m_axisMaxSliderZ->setValue(CountZ - 1);
        how_much = 2;
    }
}

void SurfaceGraph::enableResidModel(bool enable)
{
    if (enable) {
        this->nn = n;
        this->mm = m;
        m_ResidSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_ResidSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(Min, Max);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(Min, Max);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);

        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_ResidSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemeArmyBlue);

        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::red);
        m_graph->seriesList().at(0)->setBaseGradient(gr);
        m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);

        m_rangeMinX = Min;
        m_rangeMinZ = Min;
        m_stepX = (Max - Min) / double(CountX - 1);
        m_stepZ = (Max - Min) / double(CountZ - 1);
        m_axisMinSliderX->setMaximum(CountX - 2);
        m_axisMinSliderX->setValue(0);
        m_axisMaxSliderX->setMaximum(CountX - 1);
        m_axisMaxSliderX->setValue(CountX - 1);
        m_axisMinSliderZ->setMaximum(CountZ - 2);
        m_axisMinSliderZ->setValue(0);
        m_axisMaxSliderZ->setMaximum(CountZ - 1);
        m_axisMaxSliderZ->setValue(CountZ - 1);
        how_much = 3;
    }
}

void SurfaceGraph::adjustXMin(int min)
{
    double minX = m_stepX * double(min) + m_rangeMinX;

    int max = m_axisMaxSliderX->value();
    if (min >= max) {
        max = min + 1;
        m_axisMaxSliderX->setValue(max);
    }
    double maxX = m_stepX * max + m_rangeMinX;

    setAxisXRange(minX, maxX);
}

void SurfaceGraph::adjustXMax(int max)
{
    double maxX = m_stepX * double(max) + m_rangeMinX;

    int min = m_axisMinSliderX->value();
    if (max <= min) {
        min = max - 1;
        m_axisMinSliderX->setValue(min);
    }
    double minX = m_stepX * min + m_rangeMinX;

    setAxisXRange(minX, maxX);
}

void SurfaceGraph::adjustZMin(int min)
{
    double minZ = m_stepZ * double(min) + m_rangeMinZ;

    int max = m_axisMaxSliderZ->value();
    if (min >= max) {
        max = min + 1;
        m_axisMaxSliderZ->setValue(max);
    }
    double maxZ = m_stepZ * max + m_rangeMinZ;

    setAxisZRange(minZ, maxZ);
}

void SurfaceGraph::adjustZMax(int max)
{
    double maxX = m_stepZ * double(max) + m_rangeMinZ;

    int min = m_axisMinSliderZ->value();
    if (max <= min) {
        min = max - 1;
        m_axisMinSliderZ->setValue(min);
    }
    double minX = m_stepZ * min + m_rangeMinZ;

    setAxisZRange(minX, maxX);
}

void SurfaceGraph::setAxisXRange(double min, double max)
{
    m_graph->axisX()->setRange(min, max);
}

void SurfaceGraph::setAxisZRange(double min, double max)
{
    m_graph->axisZ()->setRange(min, max);
}

double function(double x, double z) {
    return qSin(x*x) + qSin(z*z);
}

double didx (double *x, double* f, int i, int j, double x0, double fx0, double xl, double fl){
    if (i == 0){
        return (f[m + j] - fx0)/(x[1] - x0);
    } else if (i == n-1){
        return (fl - f[(n-2)*m + j])/(xl - x[n-2]);
    } else {
        return (f[(i+1)*m + j] - f[(i-1)*m + j])/(x[i+1]-x[i-1]);
    }
}

double didz (double *z, double* f, int i, int j, double z0, double fz0, double zl, double fzl){
    if (j == 0){
        return (f[i*m + 1] - fz0)/(z[1] - z0);
    } else if (j == m-1){
        return (fzl - f[i*m + m-2])/(zl - z[m-2]);
    } else {
        return (f[i*m + j+1] - f[i*m + j-1])/(z[j+1] - z[j-1]);
    }
}

void transp (double *a, int u){
    double tmp;
    int i, j;
    for (i = 0; i < u-1; i++){
        for (j = i+1; j < u; j++){
            tmp = a[i*u + j];
            a[i*u + j] = a[j*u + i];
            a[j*u + i] = tmp;
        }
    }
}

void Aminus(double* A, double h) {
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 4; j++){
            if (i == j) {
                A[i*4 + j] = 1.;
            } else {
                A[i*4 + j] = 0.;
            }
        }
    }
    A[2*4] = -3./(h*h);
    A[2*4 + 1] = -2./(h);
    A[2*4 + 2] = 3./(h*h);
    A[2*4 + 3] = -1./(h);
    A[3*4] = 2./(h*h*h);
    A[3*4 + 1] = 1./(h*h);
    A[3*4 + 2] = -2./(h*h*h);
    A[3*4 + 3] = 1./(h*h);
}

void MatrMult(double* a, double* b, double* c){
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            c[i*4 + j] = a[i*4]*b[j] + a[i*4 + 1]*b[4 + j] + a[i*4 + 2]*b[2*4 + j] + a[i*4 + 3]*b[3*4 + j];
        }
    }
}

void Fij(double* f, double* dx, double* dz, double* dxdz, double* F, int i, int j){
    for (int k = 0; k < 4; k++){
        if (k % 2 == 0) {
            F[k*4] = f[(i + (int)(k/2))*m + j];
            F[k*4 + 1] = dz[(i + (int)(k/2))*m + j];
            F[k*4 + 2] = f[(i + (int)(k/2))*m + j+1];
            F[k*4 + 3] = dz[(i + (int)(k/2))*m + j+1];
        } else {
            F[k*4] = dx[(i + (int)(k/2))*m + j];
            F[k*4 + 1] = dxdz[(i + (int)(k/2))*m + j];
            F[k*4 + 2] = dx[(i + (int)(k/2))*m + j+1];
            F[k*4 + 3] = dxdz[(i + (int)(k/2))*m + j+1];
        }
    }
/*
    F[0] = f[i*m + j];
    F[1] = dz[i*m + j];
    F[2] = f[i*m + j+1];
    F[3] = dz[i*m + j+1];
    F[4] = dx[i*m + j];
    F[5] = dxdz[i*m + j];
    F[6] = dz[i*m + j+1];
    F[7] = dxdz[i*m + j+1];
    F[8] = f[(i+1)*m + j];
    F[9] = dz[(i+1)*m + j];
    F[10] = f[(i+1)*m + j+1];
    F[11] = dz[(i+1)*m + j+1];
    F[12] = dx[(i+1)*m + j];
    F[13] = dxdz[(i+1)*m + j];
    F[14] = dz[(i+1)*m + j+1];
    F[15] = dxdz[(i+1)*m + j+1];
*/}

double Pf(double**** g, double* x, double* z, int i, int j, double X, double Z){
    double y = 0.;
    for (int k = 0; k < 4; k++){
        for (int l = 0; l < 4; l++){
            y += g[i][j][k][l] * qPow(fabs(X - x[i]), k) * qPow(fabs(Z - z[j]), l);
        }
    }
    return y;
}

void SurfaceGraph::slotShortcut1() {
    how_much = (how_much + 1) % 4;
    switch (how_much){
        case 0:
            SurfaceGraph::enableInitPlotModel(true);
            break;
        case 1:
            SurfaceGraph::enableInterpModel(true);
            break;
        case 2:
            SurfaceGraph::enableInitInterpModel(true);
            break;
        case 3:
            SurfaceGraph::enableResidModel(true);
            break;
    }
    setMode(m_init, m_interp, m_together, m_resid);
}

void SurfaceGraph::slotShortcut2(){
    n *= 2;
    this->nn = n;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN, m_countM);
    switch (how_much){
        case 0:
            SurfaceGraph::enableInitPlotModel(true);
            break;
        case 1:
            SurfaceGraph::enableInterpModel(true);
            break;
        case 2:
            SurfaceGraph::enableInitInterpModel(true);
            break;
        case 3:
            SurfaceGraph::enableResidModel(true);
            break;
    }
}

void SurfaceGraph::slotShortcut4(){
    m *= 2;
    this->mm = m;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN, m_countM);
    switch (how_much){
        case 0:
            SurfaceGraph::enableInitPlotModel(true);
            break;
        case 1:
            SurfaceGraph::enableInterpModel(true);
            break;
        case 2:
            SurfaceGraph::enableInitInterpModel(true);
            break;
        case 3:
            SurfaceGraph::enableResidModel(true);
            break;
    }

}

void SurfaceGraph::slotShortcut3(){
    if ((n % 2 == 0) && (n > 4))
        n /= 2;
    this->nn = n;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN, m_countM);
    switch (how_much){
        case 0:
            SurfaceGraph::enableInitPlotModel(true);
            break;
        case 1:
            SurfaceGraph::enableInterpModel(true);
            break;
        case 2:
            SurfaceGraph::enableInitInterpModel(true);
            break;
        case 3:
            SurfaceGraph::enableResidModel(true);
            break;
    }
}

void SurfaceGraph::slotShortcut5(){
    if ((m % 2 == 0) && (m > 4))
        m /= 2;
    this->nn = n;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN, m_countM);
    switch (how_much){
        case 0:
            SurfaceGraph::enableInitPlotModel(true);
            break;
        case 1:
            SurfaceGraph::enableInterpModel(true);
            break;
        case 2:
            SurfaceGraph::enableInitInterpModel(true);
            break;
        case 3:
            SurfaceGraph::enableResidModel(true);
            break;
    }
}

void SurfaceGraph::setCount(QLabel *countN, QLabel* countM){
    countN->setNum(this->nn);
    countM->setNum(this->mm);
    m_countN = countN;
    m_countM = countM;
}

void SurfaceGraph::setMode(QRadioButton *init, QRadioButton *interp, QRadioButton *together, QRadioButton *resid){
    switch (how_much){
        case 0:
            init->setChecked(true);
            break;
        case 1:
            interp->setChecked(true);
            break;
        case 2:
            together->setChecked(true);
            break;
        case 3:
            resid->setChecked(true);
            break;
    }
    m_init = init;
    m_interp = interp;
    m_together = together;
    m_resid = resid;
}

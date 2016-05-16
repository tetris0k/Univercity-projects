#include "surfacegraph.h"

#include <QtDataVisualization/QValue3DAxis>
#include <QtDataVisualization/Q3DTheme>
#include <QtGui/QImage>
#include <QtCore/qmath.h>

#define EPS 10e-10

using namespace QtDataVisualization;

int skolkoX=40;
int skolkoZ=40;
double Min=-2.0;
double Max=2.0;
int n=4;
int m=4;
int kolichestvo=0;
double xMin=Min;
double xMax=Max;
double zMin=Min;
double zMax=Max;

SurfaceGraph::SurfaceGraph(Q3DSurface *surface)
    : m_graph(surface)
{
    m_graph->setAxisX(new QValue3DAxis);
    m_graph->setAxisY(new QValue3DAxis);
    m_graph->setAxisZ(new QValue3DAxis);
    this->nn=n;
    this->mm=m;
    m_InitPlotProxy=new QSurfaceDataProxy();
    m_InitPlotSeries=new QSurface3DSeries(m_InitPlotProxy);
    fillInitPlotProxy();
    m_InterpProxy=new QSurfaceDataProxy();
    m_InterpSeries=new QSurface3DSeries(m_InterpProxy);
    fillInterpProxy();
    m_ResidProxy=new QSurfaceDataProxy();
    m_ResidSeries=new QSurface3DSeries(m_ResidProxy);
    fillResidProxy();

}

SurfaceGraph::~SurfaceGraph()
{
    delete m_graph;
}


void SurfaceGraph::enableInitPlotModel(bool enable)
{
    if (enable){
        this->nn=n;
        this->mm=m;
        m_InitPlotSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InitPlotSeries->setFlatShadingEnabled(true);

        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(xMin, xMax);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(zMin, zMax);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);

        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_InitPlotSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemePrimaryColors);

        m_rangeMinX=xMin;
        m_rangeMinZ=zMin;
        m_stepX=(xMax-xMin)/double(skolkoX-1);
        m_stepZ=(zMax-zMin)/double(skolkoZ-1);
        kolichestvo=0;
    }
}

void SurfaceGraph::enableInterpModel(bool enable)
{
    if (enable){
        this->nn=n;
        this->mm=m;
        m_InterpSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InterpSeries->setFlatShadingEnabled(true);
        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(xMin, xMax);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(zMin, zMax);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);
        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_InterpSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemePrimaryColors);
        QLinearGradient gr;
        gr.setColorAt(0.0, Qt::green);
        m_graph->seriesList().at(0)->setBaseGradient(gr);
        m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
        m_rangeMinX=xMin;
        m_rangeMinZ=zMin;
        m_stepX=(xMax-xMin)/double(skolkoX-1);
        m_stepZ=(zMax-zMin)/double(skolkoZ-1);
        kolichestvo=1;
    }
}

void SurfaceGraph::enableInitInterpModel(bool enable)
{
    if (enable){
        this->nn=n;
        this->mm=m;
        m_InterpSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InterpSeries->setFlatShadingEnabled(true);
        m_InitPlotSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_InitPlotSeries->setFlatShadingEnabled(true);
        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(xMin,xMax);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(zMin,zMax);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);
        m_graph->setSelectionMode(QAbstract3DGraph::SelectionMultiSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_InitPlotSeries);
        m_graph->addSeries(m_InterpSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemePrimaryColors);
        QLinearGradient gr;
        gr.setColorAt(0.0,Qt::green);
        m_graph->seriesList().at(1)->setBaseGradient(gr);
        m_graph->seriesList().at(1)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
        m_rangeMinX=xMin;
        m_rangeMinZ=zMin;
        m_stepX=(xMax-xMin)/double(skolkoX-1);
        m_stepZ=(zMax-zMin)/double(skolkoZ-1);
        kolichestvo=2;
    }
}

void SurfaceGraph::enableResidModel(bool enable)
{
    if (enable){
        this->nn=n;
        this->mm=m;
        m_ResidSeries->setDrawMode(QSurface3DSeries::DrawSurfaceAndWireframe);
        m_ResidSeries->setFlatShadingEnabled(true);
        m_graph->axisX()->setLabelFormat("%.2f");
        m_graph->axisZ()->setLabelFormat("%.2f");
        m_graph->axisX()->setRange(xMin, xMax);
        m_graph->axisY()->setAutoAdjustRange(true);
        m_graph->axisZ()->setRange(zMin, zMax);
        m_graph->axisX()->setLabelAutoRotation(30);
        m_graph->axisY()->setLabelAutoRotation(90);
        m_graph->axisZ()->setLabelAutoRotation(30);
        m_graph->removeSeries(m_InterpSeries);
        m_graph->removeSeries(m_InitPlotSeries);
        m_graph->removeSeries(m_ResidSeries);
        m_graph->addSeries(m_ResidSeries);
        m_graph->activeTheme()->setType(Q3DTheme::ThemePrimaryColors);
        QLinearGradient gr;
        gr.setColorAt(0.0,Qt::green);
        m_graph->seriesList().at(0)->setBaseGradient(gr);
        m_graph->seriesList().at(0)->setColorStyle(Q3DTheme::ColorStyleRangeGradient);
        m_rangeMinX=xMin;
        m_rangeMinZ=zMin;
        m_stepX=(xMax-xMin)/double(skolkoX-1);
        m_stepZ=(zMax-zMin)/double(skolkoZ-1);
        kolichestvo=3;
    }
}

void SurfaceGraph::setAxisXRange(double min,double max)
{
    m_graph->axisX()->setRange(min,max);
}

void SurfaceGraph::setAxisZRange(double min,double max)
{
    m_graph->axisZ()->setRange(min,max);
}

void SurfaceGraph::fillInitPlotProxy()
{
    double shagX=(xMax-xMin)/double(skolkoX-1);
    double shagZ=(zMax-zMin)/double(skolkoZ-1);
    QSurfaceDataArray *dataArray=new QSurfaceDataArray;
    dataArray->reserve(skolkoZ);
    for (int i=0;i<skolkoZ;i++) {
        QSurfaceDataRow *newRow=new QSurfaceDataRow(skolkoX);
        double z=qMin(zMax,(i*shagZ+zMin));
        int index=0;
        for (int j=0;j<skolkoX;j++) {
            double x=qMin(xMax,(j*shagX+xMin));
            double y=func(x,z);
            (*newRow)[index++].setPosition(QVector3D(x,y,z));
        }
        *dataArray<<newRow;
    }
    m_InitPlotProxy->resetArray(dataArray);
}

void SurfaceGraph::fillInterpProxy()
{
    double X,Y,Z;
    int i,j,index,k,l;
    double *f=new double [n*m];
    double *x=new double [n];
    double *z=new double [m];
    double *fproizvX=new double [n*m];
    double *fproizvZ=new double [n*m];
    double *f2proizv=new double [n*m];
    double ****g=new double *** [n-1];
    for(i=0;i<n-1;i++){
        g[i]=new double ** [m-1];
        for (j=0;j<m-1;j++){
            g[i][j]=new double * [4];
            for (k=0;k<4;k++){
                g[i][j][k]=new double [4];
            }
        }
    }
    double *A=new double [16];
    double *B=new double [16];
    double *C=new double [16];
    double *F=new double [16];
    double shagXinterp=(xMax-xMin)/double(n-1);
    double shagZinterp=(zMax-zMin)/double(m-1);
    for (i=0;i<n;i++){
        x[i]=Min+i*shagXinterp;
        for (j=0;j<m;j++){
            z[j]=Min+j*shagZinterp;
            f[i*m+j]=func(x[i],z[j]);
        }
    }
    for (j=0;j<m;j++){
        for (i=0;i<n;i++){
            fproizvX[i*m+j]=proizvX(x,f,i,j);
            fproizvZ[i*m+j]=proizvZ(z,f,i,j);
        }
    }
    for (i=0;i<n;i++){
        for (j=0;j<m;j++){
             f2proizv[i*m+j]=proizvZ(z,fproizvX,i,j);
        }
    }
    for (i=0;i<n-1;i++){
        for (j=0;j<m-1;j++){
            Fij(f,fproizvX,fproizvZ, f2proizv,F,i,j);
            Aminus(A,fabs(x[i+1]-x[i]));
            Aminus(B,fabs(z[j+1]-z[j]));
            transp(B,4);
            UmnozhMatric(A,F,C);
            UmnozhMatric(C,B,F);
            for (k=0;k<4;k++){
                for (l=0;l<4;l++){
                    g[i][j][k][l]=F[k*4+l];
                }
            }
        }
    }
    int ii, jj;
    double shagX=(xMax-xMin)/double(skolkoX-1);
    double shagZ=(zMax-zMin)/double(skolkoZ-1);
    QSurfaceDataArray *dataArray=new QSurfaceDataArray;
    dataArray->reserve(skolkoZ);
    for (j=0;j<skolkoZ;j++) {
        QSurfaceDataRow *newRow=new QSurfaceDataRow(skolkoX);
        Z=qMin(zMax,(j*shagZ+zMin));
        jj=gde(Z,z,m);
        index=0;
        for (i=0;i<skolkoX;i++){
            X=qMin(xMax,(i*shagX+xMin));
            ii=gde(X,x,n);
            Y=Pf(g,x,z,ii,jj,X,Z);
            (*newRow)[index++].setPosition(QVector3D(X,Y,Z));
        }
        *dataArray<<newRow;
    }
    for(i=0;i<n-1;i++){
        for(j=0;j<m-1;j++){
            for(k=0;k<4;k++){
                delete [] g[i][j][k];
            }
            delete [] g[i][j];
        }
        delete [] g[i];
    }
    delete []g;
    delete []f;
    delete []fproizvX;
    delete []fproizvZ;
    delete []f2proizv;
    delete []A;
    delete []B;
    delete []C;
    delete []F;
    delete []x;
    delete []z;
    m_InterpProxy->resetArray(dataArray);
}

void SurfaceGraph::fillResidProxy(){
    this->nn=n;
    this->mm=m;
    double X,Y,Z;
    int i,j,index,k,l;
    double *f=new double [n*m];
    double *x=new double [n];
    double *z=new double [m];
    double *fproizvX=new double [n*m];
    double *fproizvZ=new double [n*m];
    double *f2proizv=new double [n*m];
    double ****g=new double *** [n-1];
    for(i=0;i<n-1;i++){
        g[i]=new double ** [m-1];
        for (j=0;j<m-1;j++){
            g[i][j]=new double * [4];
            for(k=0;k<4;k++){
                g[i][j][k]=new double [4];
            }
        }
    }
    double *A=new double [16];
    double *B=new double [16];
    double *C=new double [16];
    double *F=new double [16];
    double shagXinterp=(xMax-xMin)/double(n-1);
    double shagZinterp=(zMax-zMin)/double(m-1);
    for (i=0;i<n;i++){
        x[i]=Min+i*shagXinterp;
        for(j=0;j<m;j++){
            z[j]=Min+j*shagZinterp;
            f[i*m+j]=func(x[i],z[j]);
        }
    }
    for(j=0;j<m;j++){
        for(i=0;i<n;i++){
           fproizvX[i*m+j]=proizvX(x,f,i,j);
           fproizvZ[i*m+j]=proizvZ(z,f,i,j);
        }
    }
    for(i=0;i<n;i++){
        for(j=0;j<m;j++){
            f2proizv[i*m+j]=proizvZ(z,fproizvX,i,j);
        }
    }
    for(i=0;i<n-1;i++){
        for(j=0;j<m-1;j++){
            Fij(f,fproizvX,fproizvZ,f2proizv,F,i,j);
            Aminus(A,fabs(x[i+1]-x[i]));
            Aminus(B,fabs(z[j+1]-z[j]));
            transp(B,4);
            UmnozhMatric(A,F,C);
            UmnozhMatric(C,B,F);
            for(k=0;k<4;k++){
                for(l=0;l<4;l++){
                    g[i][j][k][l]=F[k*4+l];
                }
            }
        }
    }
    int ii, jj;
    double shagX=(xMax-xMin)/double(skolkoX-1);
    double shagZ=(zMax-zMin)/double(skolkoZ-1);
    QSurfaceDataArray *dataArray=new QSurfaceDataArray;
    dataArray->reserve(skolkoZ);
    for(j=0;j<skolkoZ;j++){
        QSurfaceDataRow *newRow=new QSurfaceDataRow(skolkoX);
        Z=qMin(zMax,(j*shagZ+zMin));
        jj=gde(Z,z,m);
        index=0;
        for(i=0;i<skolkoX;i++){
            X=qMin(xMax,(i*shagX+xMin));
            ii=gde(X,x,n);
            Y=fabs(Pf(g,x,z,ii,jj,X,Z)-func(X,Z));
            (*newRow)[index++].setPosition(QVector3D(X,Y,Z));
        }
        *dataArray<<newRow;
    }
    delete [] x;
    delete [] z;
    delete [] f;
    delete [] fproizvX;
    delete [] fproizvZ;
    delete [] f2proizv;
    delete [] A;
    delete [] B;
    delete [] C;
    delete [] F;
    for(i=0;i<n-1;i++){
        for(j=0;j<m-1;j++){
            for(k=0;k<4;k++){
                delete [] g[i][j][k];
            }
            delete [] g[i][j];
        }
        delete [] g[i];
    }
    delete [] g;
    m_ResidProxy->resetArray(dataArray);
}

int gde(double X,double* x,int l){
    if(X<x[0]){
        return 0;
    }
    for(int k=0;k<l-1;k++){
        if((X-x[k]>-EPS)&&(x[k+1]-X>-EPS)){
            return k;
        }
    }
    return l-2;
}

double func(double x,double z){
    return qSin(x*x)+qSin(z*z);
}

double proizvX(double *x,double* f,int i,int j){
    if(i==0){
        return (3.*(f[m+j]-f[j])/(x[1]-x[0])-proizvX(x,f,1,j))/2.;
    } else if(i==n-1){
        return (3.*(f[(n-1)*m+j]-f[(n-2)*m+j])/(x[n-1]-x[n-2])-proizvX(x,f,n-2,j))/2.;
    } else {
        return (f[(i+1)*m+j]-f[(i-1)*m+j])/(x[i+1]-x[i-1]);
    }
}


double proizvZ(double *z,double* f,int i,int j){
    if(j==0){
        return (3.*(f[i*m+1]-f[i*m])/(z[1]-z[0])+proizvZ(z,f,i,1))/2.;
    } else if(j==m-1){
        return (3.*(f[i*m+m-1]-f[i*m+m-2])/(z[m-1]-z[m-2])-proizvZ(z,f,i,m-2))/2.;
    } else {
        return (f[i*m+j+1]-f[i*m+j-1])/(z[j+1]-z[j-1]);
    }
}


void transp(double *a,int u){
    double tmp;
    int i,j;
    for(i=0;i<u-1;i++){
        for(j=i+1;j<u;j++){
            tmp=a[i*u+j];
            a[i*u+j]=a[j*u+i];
            a[j*u+i]=tmp;
        }
    }
}

void Aminus(double* A,double h) {
    for(int i=0;i<2;i++){
        for(int j=0;j<4;j++){
            if(i==j){
                A[i*4+j]=1.;
            } else {
                A[i*4+j]=0.;
            }
        }
    }
    A[2*4]=-3./(h*h);
    A[2*4+1]=-2./(h);
    A[2*4+2]=3./(h*h);
    A[2*4+3]=-1./(h);
    A[3*4]=2./(h*h*h);
    A[3*4+1]=1./(h*h);
    A[3*4+2]=-2./(h*h*h);
    A[3*4+3]=1./(h*h);
}

void UmnozhMatric(double* a,double* b,double* c){
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            c[i*4+j]=a[i*4]*b[j]+a[i*4+1]*b[4+j]+a[i*4+2]*b[2*4+j]+a[i*4+3]*b[3*4+j];
        }
    }
}

void Fij(double* f,double* dx,double* dz,double* dxdz,double* F,int i,int j){
    for(int k=0;k<4;k++){
        if(k%2==0){
            F[k*4]=f[(i+(int)(k/2))*m+j];
            F[k*4+1]=dz[(i+(int)(k/2))*m+j];
            F[k*4+2]=f[(i+(int)(k/2))*m+j+1];
            F[k*4+3]=dz[(i+(int)(k/2))*m+j+1];
        } else {
            F[k*4]=dx[(i+(int)(k/2))*m+j];
            F[k*4+1]=dxdz[(i+(int)(k/2))*m+j];
            F[k*4+2]=dx[(i+(int)(k/2))*m+j+1];
            F[k*4+3]=dxdz[(i+(int)(k/2))*m+j+1];
        }
    }
}

double Pf(double**** g,double* x,double* z,int i,int j,double X,double Z){
    double y=0.;
    for(int k=0;k<4;k++){
        for (int l=0;l<4;l++){
            y=y+g[i][j][k][l]*qPow(fabs(X-x[i]),k)*qPow(fabs(Z-z[j]),l);
        }
    }
    return y;
}

void SurfaceGraph::slotShortcut1() {
    kolichestvo=(kolichestvo+1)%4;
    setKol(m_first,m_second,m_third,m_fourth);             //при нажатии 1 используем функцию выставления точек
    if (kolichestvo==0){
        SurfaceGraph::enableInitPlotModel(true);
    } else if(kolichestvo==1){
        SurfaceGraph::enableInterpModel(true);
    } else if(kolichestvo==2){
        SurfaceGraph::enableInitInterpModel(true);
    } else {
        SurfaceGraph::enableResidModel(true);
    }
}

void SurfaceGraph::slotShortcut2(){
    n*=2;
    this->nn=n;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN,m_countM);
    if (kolichestvo==0){
        SurfaceGraph::enableInitPlotModel(true);
    } else if(kolichestvo==1){
        SurfaceGraph::enableInterpModel(true);
    } else if(kolichestvo==2){
        SurfaceGraph::enableInitInterpModel(true);
    } else {
        SurfaceGraph::enableResidModel(true);
    }
}

void SurfaceGraph::slotShortcut3(){
    m*=2;
    this->mm=m;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN,m_countM);
    if (kolichestvo==0){
        SurfaceGraph::enableInitPlotModel(true);
    } else if(kolichestvo==1){
        SurfaceGraph::enableInterpModel(true);
    } else if(kolichestvo==2){
        SurfaceGraph::enableInitInterpModel(true);
    } else {
        SurfaceGraph::enableResidModel(true);
    }
}

void SurfaceGraph::slotShortcut4(){
    xMax=xMax+0.1;
    xMin=xMin-0.1;
    m_graph->axisX()->setRange(xMin,xMax);
}

void SurfaceGraph::slotShortcut5(){
    if (xMax-0.1>xMin+0.11){
        xMax=xMax-0.1;
        xMin=xMin+0.1;
    }
    m_graph->axisX()->setRange(xMin,xMax);
}

void SurfaceGraph::slotShortcut6(){
    zMax=zMax+0.1;
    zMin=zMin-0.1;
    m_graph->axisZ()->setRange(zMin,zMax);
}

void SurfaceGraph::slotShortcut7(){
    if (zMax-0.1>zMin+0.11){
        zMax=zMax-0.1;
        zMin=zMin+0.1;
    }
    m_graph->axisZ()->setRange(zMin,zMax);
}

void SurfaceGraph::setCount(QLabel *countN,QLabel *countM){
    countN->setNum(this->nn);
    countM->setNum(this->mm);
    m_countN=countN;
    m_countM=countM;
}

void SurfaceGraph::setKol(QRadioButton *first, QRadioButton *second,QRadioButton *third,QRadioButton *fourth){
    if (kolichestvo == 0){
        first->setChecked(true);
    } else if (kolichestvo == 1){
        second->setChecked(true);
    } else if (kolichestvo==2){
        third->setChecked(true);    //Эта функция делает отмеченными кнопочки. работает, как setNum (чуть выше)
    } else {
        fourth->setChecked(true);
    }
    m_first=first;                  //выставляем нужные точки, и синхронизируем их с добавленными в класс SurfaceGraph (как в setNum)
    m_second=second;
    m_third=third;
    m_fourth=fourth;
}

void SurfaceGraph::slotShortcut8(){
    if ((n%2 == 0)&&(n > 2))
        n/=2;
    this->nn=n;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN,m_countM);
    if (kolichestvo==0){
        SurfaceGraph::enableInitPlotModel(true);
    } else if(kolichestvo==1){
        SurfaceGraph::enableInterpModel(true);
    } else if(kolichestvo==2){
        SurfaceGraph::enableInitInterpModel(true);
    } else {
        SurfaceGraph::enableResidModel(true);
    }
}

void SurfaceGraph::slotShortcut9(){
    if ((m%2==0)&&(m>2))
        m/=2;
    this->mm=m;
    fillInitPlotProxy();
    fillInterpProxy();
    fillResidProxy();
    setCount(m_countN,m_countM);
    if (kolichestvo==0){
        SurfaceGraph::enableInitPlotModel(true);
    } else if(kolichestvo==1){
        SurfaceGraph::enableInterpModel(true);
    } else if(kolichestvo==2){
        SurfaceGraph::enableInitInterpModel(true);
    } else {
        SurfaceGraph::enableResidModel(true);
    }
}

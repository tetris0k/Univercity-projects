#include "surfacegraph.h"

#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLabel>
#include <QtGui/QPainter>
#include <QtGui/QScreen>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    Q3DSurface *graph = new Q3DSurface();
    QWidget *container = QWidget::createWindowContainer(graph);

    QSize screenSize = graph->screen()->size();
    container->setMinimumSize(QSize(screenSize.width() / 2, screenSize.height() / 1.6));
    container->setMaximumSize(screenSize);
    container->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    container->setFocusPolicy(Qt::StrongFocus);

    QWidget *widget = new QWidget;
    QHBoxLayout *hLayout = new QHBoxLayout(widget);
    QVBoxLayout *vLayout = new QVBoxLayout();
    hLayout->addWidget(container, 1);
    hLayout->addLayout(vLayout);
    vLayout->setAlignment(Qt::AlignTop);

    widget->setWindowTitle(QStringLiteral("Interpolation"));

    QGroupBox *modelGroupBox = new QGroupBox(QStringLiteral("Mode"));

    QRadioButton *InitPlotModelRB = new QRadioButton(widget);
    InitPlotModelRB->setText(QStringLiteral("Initial plot"));
    InitPlotModelRB->setChecked(false);

    QRadioButton *InterpModelRB = new QRadioButton(widget);
    InterpModelRB->setText(QStringLiteral("Interpolating plot"));
    InterpModelRB->setChecked(false);

    QRadioButton *InitInterpModelRB = new QRadioButton(widget);
    InitInterpModelRB->setText(QStringLiteral("Initial and Interpolating plots"));
    InitInterpModelRB->setChecked(false);

    QRadioButton *ResidModelRB = new QRadioButton(widget);
    ResidModelRB->setText(QStringLiteral("Residual norm plot"));
    ResidModelRB->setChecked(false);

    QVBoxLayout *modelVBox = new QVBoxLayout;
    modelVBox->addWidget(InitPlotModelRB);
    modelVBox->addWidget(InterpModelRB);
    modelVBox->addWidget(InitInterpModelRB);
    modelVBox->addWidget(ResidModelRB);
    modelGroupBox->setLayout(modelVBox);

    QSlider *axisMinSliderX = new QSlider(Qt::Horizontal, widget);
    axisMinSliderX->setMinimum(0);
    axisMinSliderX->setTickInterval(1);
    axisMinSliderX->setEnabled(true);
    QSlider *axisMaxSliderX = new QSlider(Qt::Horizontal, widget);
    axisMaxSliderX->setMinimum(1);
    axisMaxSliderX->setTickInterval(1);
    axisMaxSliderX->setEnabled(true);
    QSlider *axisMinSliderZ = new QSlider(Qt::Horizontal, widget);
    axisMinSliderZ->setMinimum(0);
    axisMinSliderZ->setTickInterval(1);
    axisMinSliderZ->setEnabled(true);
    QSlider *axisMaxSliderZ = new QSlider(Qt::Horizontal, widget);
    axisMaxSliderZ->setMinimum(1);
    axisMaxSliderZ->setTickInterval(1);
    axisMaxSliderZ->setEnabled(true);

    QLabel *countN = new QLabel(widget);
    QLabel *countM = new QLabel(widget);

    vLayout->addWidget(modelGroupBox);
    vLayout->addWidget(new QLabel(QStringLiteral("Column range")));
    vLayout->addWidget(axisMinSliderX);
    vLayout->addWidget(axisMaxSliderX);
    vLayout->addWidget(new QLabel(QStringLiteral("Row range")));
    vLayout->addWidget(axisMinSliderZ);
    vLayout->addWidget(axisMaxSliderZ);
    vLayout->addWidget(new QLabel(QStringLiteral("Points at X axis")));
    vLayout->addWidget(countN);
    vLayout->addWidget(new QLabel(QStringLiteral("Points at Z axis")));
    vLayout->addWidget(countM);
    vLayout->addWidget(new QLabel(QStringLiteral("Press 1 to change the plot\n"
                                                 "Press 2 to increase number of points at X axis\n"
                                                 "Press 3 to decrease number of points at X axis\n"
                                                 "Press 4 to increase number of points at Z axis\n"
                                                 "Press 5 to decrease number of points at Z axis")));

    widget->show();

    SurfaceGraph *modifier = new SurfaceGraph(graph);

    modifier->key1 = new QShortcut(widget);
    modifier->key1->setKey(Qt::Key_1);
    QObject::connect(modifier->key1, SIGNAL(activated()), modifier, SLOT(slotShortcut1()));

    modifier->key2 = new QShortcut(widget);
    modifier->key2->setKey(Qt::Key_2);
    QObject::connect(modifier->key2, SIGNAL(activated()), modifier, SLOT(slotShortcut2()));

    modifier->key3 = new QShortcut(widget);
    modifier->key3->setKey(Qt::Key_3);
    QObject::connect(modifier->key3, SIGNAL(activated()), modifier, SLOT(slotShortcut3()));

    modifier->key4 = new QShortcut(widget);
    modifier->key4->setKey(Qt::Key_4);
    QObject::connect(modifier->key4, SIGNAL(activated()), modifier, SLOT(slotShortcut4()));

    modifier->key5 = new QShortcut(widget);
    modifier->key5->setKey(Qt::Key_5);
    QObject::connect(modifier->key5, SIGNAL(activated()), modifier, SLOT(slotShortcut5()));

    QObject::connect(InitPlotModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableInitPlotModel);
    QObject::connect(InterpModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableInterpModel);
    QObject::connect(InitInterpModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableInitInterpModel);
    QObject::connect(ResidModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableResidModel);
    QObject::connect(axisMinSliderX, &QSlider::valueChanged,
                     modifier, &SurfaceGraph::adjustXMin);
    QObject::connect(axisMaxSliderX, &QSlider::valueChanged,
                     modifier, &SurfaceGraph::adjustXMax);
    QObject::connect(axisMinSliderZ, &QSlider::valueChanged,
                     modifier, &SurfaceGraph::adjustZMin);
    QObject::connect(axisMaxSliderZ, &QSlider::valueChanged,
                     modifier, &SurfaceGraph::adjustZMax);

    modifier->setAxisMinSliderX(axisMinSliderX);
    modifier->setAxisMaxSliderX(axisMaxSliderX);
    modifier->setAxisMinSliderZ(axisMinSliderZ);
    modifier->setAxisMaxSliderZ(axisMaxSliderZ);

    InitPlotModelRB->setChecked(true);
    modifier->setCount(countN, countM);
    modifier->setMode(InitPlotModelRB, InterpModelRB, InitInterpModelRB, ResidModelRB);

    return app.exec();
}

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
    hLayout->addLayout(vLayout);
    hLayout->addWidget(container, 1);
    vLayout->setAlignment(Qt::AlignTop);
    widget->setWindowTitle(QStringLiteral("Интерполяция"));
    QGroupBox *modelGroupBox = new QGroupBox(QStringLiteral("Графики"));
    QRadioButton *InitPlotModelRB = new QRadioButton(widget);
    InitPlotModelRB->setText(QStringLiteral("Начальный график"));
    InitPlotModelRB->setChecked(false);
    QRadioButton *InterpModelRB = new QRadioButton(widget);
    InterpModelRB->setText(QStringLiteral("Интерполяция"));
    InterpModelRB->setChecked(false);
    QRadioButton *InitInterpModelRB = new QRadioButton(widget);
    InitInterpModelRB->setText(QStringLiteral("Изначальный и интерполяция"));
    InitInterpModelRB->setChecked(false);
    QRadioButton *ResidModelRB = new QRadioButton(widget);
    ResidModelRB->setText(QStringLiteral("Невязка"));
    ResidModelRB->setChecked(false);
    QVBoxLayout *modelVBox = new QVBoxLayout;
    modelVBox->addWidget(InitPlotModelRB);
    modelVBox->addWidget(InterpModelRB);
    modelVBox->addWidget(InitInterpModelRB);
    modelVBox->addWidget(ResidModelRB);
    modelGroupBox->setLayout(modelVBox);
    QLabel *countN = new QLabel(widget);
    QLabel *countM = new QLabel(widget);
    vLayout->addWidget(modelGroupBox);
    vLayout->addWidget(new QLabel(QStringLiteral("точек на оси Х")));
    vLayout->addWidget(countN);
    vLayout->addWidget(new QLabel(QStringLiteral("точек на оси Z")));
    vLayout->addWidget(countM);
    vLayout->addWidget((new QLabel(QStringLiteral("1-графики\n"
                                                  "2-больше точек на Х\n"
                                                  "3-больше точек на Z\n"
                                                  "4-увеличить диапазон X\n"
                                                  "5-уменьшить диапазон X\n"
                                                  "6-увеличить диапазон Z\n"
                                                  "7-уменьшить диапазон Z\n"
                                                  "8-меньше точек на X\n"
                                                  "9-меньше точек на Z"))));
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
    modifier->key6 = new QShortcut(widget);
    modifier->key6->setKey(Qt::Key_6);
    QObject::connect(modifier->key6, SIGNAL(activated()), modifier, SLOT(slotShortcut6()));
    modifier->key7 = new QShortcut(widget);
    modifier->key7->setKey(Qt::Key_7);
    QObject::connect(modifier->key7, SIGNAL(activated()), modifier, SLOT(slotShortcut7()));
    modifier->key8 = new QShortcut(widget);
    modifier->key8->setKey(Qt::Key_8);
    QObject::connect(modifier->key8, SIGNAL(activated()), modifier, SLOT(slotShortcut8()));
    modifier->key9 = new QShortcut(widget);
    modifier->key9->setKey(Qt::Key_9);
    QObject::connect(modifier->key9, SIGNAL(activated()), modifier, SLOT(slotShortcut9()));


    QObject::connect(InitPlotModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableInitPlotModel);
    QObject::connect(InterpModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableInterpModel);
    QObject::connect(InitInterpModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableInitInterpModel);
    QObject::connect(ResidModelRB, &QRadioButton::toggled,
                     modifier, &SurfaceGraph::enableResidModel);
    InitPlotModelRB->setChecked(true);
    modifier->setCount(countN, countM);
    modifier->setKol(InitPlotModelRB,InterpModelRB, InitInterpModelRB,ResidModelRB);  //изначально выставим точки
    return app.exec();
}

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QVector3D>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits{HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits>  MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    double contangente(double angle);
    float AireBarycentrique(MyMesh* _mesh, int vertexID);

    float faceArea(MyMesh* _mesh, int faceID);
    QVector3D LaplaceBeltrami(MyMesh* _mesh, int vertexID);
    MyMesh meshLoaded;
    int VertexLaplace = 0;


private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

    void on_vertexSelec_valueChanged(const QString &arg1);

    void on_laplace_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H

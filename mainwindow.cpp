#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    // Cet exemple montre comment afficher un cube (sans passer par une structure de maillage) directement dans un MeshViewerWidget

    // une liste de sommets (x, y, z)
    GLfloat verts[24] = {
        -1.0f, 1.0f, -1.0f,
        -1.0f, -1.0f, -1.0f,
        -1.0f, 1.0f, 1.0f,
        -1.0f, -1.0f, 1.0f,
        1.0f, 1.0f, 1.0f,
        1.0f, -1.0f, 1.0f,
        1.0f, 1.0f, -1.0f,
        1.0f, -1.0f, -1.0f
    };

    // une liste de couleur (r, g, b)
    GLfloat cols[24] = {
        0.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f,
        1.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 1.0f
    };

    // une liste de triangles coresspondants aux sommets de verts
    GLuint IndiceArray[36] = {
        0,1,2,    2,1,3,
        4,5,6,    6,5,7,
        3,1,5,    5,1,7,
        0,2,6,    6,2,4,
        6,7,0,    0,7,1,
        2,3,4,    4,3,5
    };

    ui->widget->loadMesh(verts, cols, 24, IndiceArray, 36);
}

void MainWindow::on_pushButton_2_clicked()
{
    // Cet exemple montre comment ouvrir un .obj, et l'afficher dans le MeshViewerWidget
    MyMesh mesh;
    OpenMesh::IO::read_mesh(mesh, "../meshFiles/bunnyLowPoly.obj");

    GLuint* IndiceArray = new GLuint[mesh.n_faces() * 3];

    MyMesh::ConstFaceIter fIt(mesh.faces_begin()), fEnd(mesh.faces_end());
    MyMesh::ConstFaceVertexIter fvIt;
    int i = 0;
    for (; fIt!=fEnd; ++fIt)
    {
        fvIt = mesh.cfv_iter(*fIt);
        IndiceArray[i] = fvIt->idx(); i++; ++fvIt;
        IndiceArray[i] = fvIt->idx(); i++; ++fvIt;
        IndiceArray[i] = fvIt->idx(); i++;
    }

    GLfloat* cols = new GLfloat[mesh.n_vertices() * 3];
    for(int c = 0; c < mesh.n_vertices() * 3; c = c + 3)
    {
        cols[c] = 0.5; cols[c+1] = 0.5; cols[c+2] = 0.5;
    }

    ui->widget->loadMesh((GLfloat*)mesh.points(), cols, mesh.n_vertices() * 3, IndiceArray, mesh.n_faces() * 3);
}

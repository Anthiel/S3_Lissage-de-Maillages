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
    meshLoaded = mesh;
    ui->widget->loadMesh((GLfloat*)mesh.points(), cols, mesh.n_vertices() * 3, IndiceArray, mesh.n_faces() * 3);
}

float MainWindow::AireBarycentrique(MyMesh* _mesh, int vertexID){
    float aireTotal;
     VertexHandle v_it = _mesh->vertex_handle(vertexID);

     // parcours des faces autour de vertexID
    for(MyMesh::VertexFaceIter  vf_it = _mesh->vf_iter(v_it); vf_it; ++vf_it) {
        FaceHandle fh = *vf_it;
        aireTotal += faceArea(_mesh, fh.idx());
    }
    return 1/3.0*aireTotal;
}


float MainWindow::faceArea(MyMesh* _mesh, int faceID)
{
    FaceHandle fh = _mesh->face_handle(faceID);
    std::vector <int> pointID;
    VectorT <float,3> points[3];

    for (MyMesh::FaceVertexIter curVert = _mesh->fv_iter(fh); curVert.is_valid(); curVert ++)
    {
        pointID.push_back((*curVert).idx());
    }

    for (int i=0; i<pointID.size();i++)
    {
        VertexHandle vh = _mesh->vertex_handle(pointID.at(i));
        points[i][0] = _mesh->point(vh)[0];
        points[i][1] = _mesh->point(vh)[1];
        points[i][2] = _mesh->point(vh)[2];
    }
    VectorT<float,3> BmoinsA = points[1] - points[0];
    VectorT<float,3> CmoinsA = points[2] - points[0];

    VectorT<float,3> res;
    res[0] = BmoinsA[1] * CmoinsA[2] - BmoinsA[2] * CmoinsA[1];
    res[1] = BmoinsA[2] * CmoinsA[0] - BmoinsA[0] * CmoinsA[2];
    res[2] = BmoinsA[0] * CmoinsA[1] - BmoinsA[1] * CmoinsA[0];

    return res.norm()/2.0;
}

QVector3D MainWindow::LaplaceBeltrami(MyMesh* _mesh, int vertexID){

    VertexHandle v_it = _mesh->vertex_handle(vertexID);
    std::vector<QVector3D> AngleID;

    MyMesh::Point sum;

    //parcours des vertexs autour du vertex
    for(MyMesh::VertexVertexIter  vv_it = _mesh->vv_iter(v_it); vv_it; ++vv_it) {
        VertexHandle vh = *vv_it;
        AngleID.clear();

        //on cherche les deux faces des deux vertexs
        for (MyMesh::VertexFaceIter curVert = _mesh->vf_iter(v_it); curVert.is_valid(); curVert ++)
        {
           for (MyMesh::VertexFaceIter curVert1 = _mesh->vf_iter(vh); curVert1.is_valid(); curVert1 ++)
           {
               if((*curVert).idx() == (*curVert1).idx()){
               //cette face a les deux vertexs
               for(MyMesh::FaceVertexIter fh = _mesh->fv_iter(curVert); fh.is_valid(); fh ++){
                   //on parcours les vertexs de cette face pour trouver le vertex manquant
                   VertexHandle vh1 = *fh;
                   if(vh1.idx() != vh.idx() && vh1.idx() != v_it.idx()){
                       AngleID.push_back(QVector3D(vh.idx(), vh1.idx() , v_it.idx()));
                   }
               }
            }
         }
     }
     //arcos(produit scalaire entre deux vecteurs) = angle entre deux vecteurs

     //premier angle
     MyMesh::Point pointV = _mesh->point(_mesh->vertex_handle(vertexID));
     MyMesh::Point pointVi = _mesh->point(_mesh->vertex_handle(vh.idx()));

     MyMesh::Point pointCentral = _mesh->point(_mesh->vertex_handle(AngleID.at(0).y()));
     MyMesh::Point pointA = _mesh->point(_mesh->vertex_handle(AngleID.at(0).x()));
     MyMesh::Point pointB = _mesh->point(_mesh->vertex_handle(AngleID.at(0).z()));

     MyMesh::Point vecteurCentralA = pointCentral - pointA;
     MyMesh::Point vecteurCentralB = pointCentral - pointB;

     double angleAlphaI = acos(vecteurCentralA | vecteurCentralB);

     pointCentral = _mesh->point(_mesh->vertex_handle(AngleID.at(1).y()));
     pointA = _mesh->point(_mesh->vertex_handle(AngleID.at(1).x()));
     pointB = _mesh->point(_mesh->vertex_handle(AngleID.at(1).z()));

     vecteurCentralA = pointCentral - pointA;
     vecteurCentralB = pointCentral - pointB;

     double angleBetaI = acos(vecteurCentralA | vecteurCentralB);

     sum += ( contangente(angleAlphaI) + contangente(angleBetaI) ) * (pointVi - pointV);
    }

    MyMesh::Point res = 1/(2*AireBarycentrique(_mesh, vertexID)) * sum;
    _mesh->set_point(v_it, res);
}

double MainWindow::contangente(double angle){
    return 1/tan(angle);
}

void MainWindow::on_vertexSelec_valueChanged(const QString &arg1)
{
    VertexLaplace = arg1.toInt();
}

void MainWindow::on_laplace_clicked()
{
    LaplaceBeltrami(&meshLoaded, VertexLaplace);
}

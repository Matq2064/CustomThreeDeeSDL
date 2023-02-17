#define SDL_MAIN_HANDLED
#include <SDL.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

#include "Clock.h"
#include "Colors.h"

template <class T>
T clampLower(T x, T minimum) {
    if (x < minimum)
        return minimum;

    return x;
}

struct vec3 {
    double x, y, z;
};

struct vec2 {
    double x, y;
};

struct TrianglePlane {
    vec3 points[3];
};

struct Triangle {
    vec2 m_aPoints[3];

    Triangle() : m_aPoints() {}
    explicit Triangle(TrianglePlane Projected) : m_aPoints() {
        m_aPoints[0].x = Projected.points[0].x;
        m_aPoints[0].y = Projected.points[0].y;
        m_aPoints[1].x = Projected.points[1].x;
        m_aPoints[1].y = Projected.points[1].y;
        m_aPoints[2].x = Projected.points[2].x;
        m_aPoints[2].y = Projected.points[2].y;
    }
};

struct mesh {
    std::vector<TrianglePlane> triangles;
};

struct Matrix {
    double values[4][4] = {0};
};

void SetProjectionMatrix(Matrix& MatrixOut, double Near, double Far, double FieldOfView, double AspectRatio) {
    double FieldOfViewRadians = 1.0 / tan(FieldOfView * 0.5 / 180.0 * M_PI);

    MatrixOut.values[0][0] = AspectRatio * FieldOfViewRadians;
    MatrixOut.values[1][1] = FieldOfViewRadians;
    MatrixOut.values[2][2] = Far / (Far - Near);
    MatrixOut.values[3][2] = (-Far * Near) / (Far - Near);
    MatrixOut.values[2][3] = 1.0;
    MatrixOut.values[3][3] = 0.0;
}

void SetRotationXMatrix(Matrix& MatrixOut, double radians) {
    double cos = std::cos(radians);
    double sin = std::sin(radians);
    MatrixOut.values[0][0] = 1.0;
    MatrixOut.values[1][1] = cos;
    MatrixOut.values[1][2] = sin;
    MatrixOut.values[2][1] = -sin;
    MatrixOut.values[2][2] = cos;
    MatrixOut.values[3][3] = 1.0;
}

void SetRotationYMatrix(Matrix& MatrixOut, double radians) {
    double cos = std::cos(radians);
    double sin = std::sin(radians);
    MatrixOut.values[0][0] = cos;
    MatrixOut.values[0][2] = -sin;
    MatrixOut.values[1][1] = 1.0;
    MatrixOut.values[2][0] = sin;
    MatrixOut.values[2][2] = cos;
    MatrixOut.values[3][3] = 1.0;
}

void SetRotationZMatrix(Matrix& MatrixOut, double radians) {
    double cos = std::cos(radians);
    double sin = std::sin(radians);
    MatrixOut.values[0][0] = cos;
    MatrixOut.values[0][1] = sin;
    MatrixOut.values[1][0] = -sin;
    MatrixOut.values[1][1] = cos;
    MatrixOut.values[2][2] = 1.0;
    MatrixOut.values[3][3] = 1.0;
}

void PointMove(vec2& Point, double x, double y) {
    Point.x += x;
    Point.y += y;
}

void PointScale(vec2& Point, double mul_x, double mul_y) {
    Point.x *= mul_x;
    Point.y *= mul_y;
}

void PointApplyMatrix(vec3& PointIn, vec3& PointOut, Matrix& matrix) {
    PointOut.x = PointIn.x * matrix.values[0][0] + PointIn.y * matrix.values[1][0] + PointIn.z * matrix.values[2][0] + matrix.values[3][0];
    PointOut.y = PointIn.x * matrix.values[0][1] + PointIn.y * matrix.values[1][1] + PointIn.z * matrix.values[2][1] + matrix.values[3][1];
    PointOut.z = PointIn.x * matrix.values[0][2] + PointIn.y * matrix.values[1][2] + PointIn.z * matrix.values[2][2] + matrix.values[3][2];
    double w = PointIn.x * matrix.values[0][3] + PointIn.y * matrix.values[1][3] + PointIn.z * matrix.values[2][3] + matrix.values[3][3];

    if (w != 0.0)
        PointOut.x /= w; PointOut.y /= w; PointOut.z /= w;
}

// void DrawLine(SDL_Renderer* renderer, double x1, double y1, double x2, double y2) {
//     SDL_RenderDrawLine(renderer, (int)(x1), (int)(y1), (int)(x2), (int)(y2));
// }

double GetDistance(vec2 Vector) {
    return std::sqrt(std::pow(Vector.x, 2) + std::pow(Vector.y, 2));
}

double GetDistance(vec2 Point1, vec2 Point2) {
    return std::sqrt(std::pow(Point2.x - Point1.x, 2) + std::pow(Point2.y - Point1.y, 2));
}

void DrawLine(SDL_Renderer* renderer, vec2 Point1, vec2 Point2) {
    vec2 Travel = { Point2.x - Point1.x, Point2.y - Point1.y };
    double Distance = GetDistance(Travel);
    vec2 Slice = { Travel.x / Distance, Travel.y / Distance };

    int Iterations = int(Distance + 1);
    for (int i = 0; i < Iterations; i++) {
        vec2 CurrentPoint = { Point1.x + Slice.x * i, Point1.y + Slice.y * i };
        SDL_RenderDrawPoint(renderer, int(std::round(CurrentPoint.x)), int(std::round(CurrentPoint.y)));
    }
    SDL_RenderDrawPoint(renderer, int(std::round(Point2.x)), int(std::round(Point2.y)));
}

void DrawTriangle(SDL_Renderer* renderer, Triangle& triangle) {
    DrawLine(renderer, triangle.m_aPoints[0], triangle.m_aPoints[1]);
    DrawLine(renderer, triangle.m_aPoints[1], triangle.m_aPoints[2]);
    DrawLine(renderer, triangle.m_aPoints[2], triangle.m_aPoints[0]);
}

std::vector<int> Interpolate(const vec2& Point1, const vec2& Point2) {
    std::vector<int> xCoords;
    // first and last coordinate
    int Height = (int)(Point2.y) - (int)(Point1.y);
    double TotalX = (Point2.x - Point1.x);
    double DeltaX = TotalX / (Point2.y - Point1.y);

    double CurrentX = Point1.x;
    xCoords.push_back((int)(round(CurrentX))); // First point

    for (int i = 0; i < Height-1; i++)
        xCoords.push_back((int)(round(CurrentX + (i + 1) * DeltaX)));

    xCoords.push_back((int)(round(CurrentX + TotalX))); // Last point

    return xCoords;
}

void FillTriangle(SDL_Renderer* Renderer, Triangle& triangle) {
    vec2 Point0 = triangle.m_aPoints[0];
    vec2 Point1 = triangle.m_aPoints[1];
    vec2 Point2 = triangle.m_aPoints[2];

    // Sort the m_aPoints so that y0 <= y1 <= y2
    if (Point1.y < Point0.y) { std::swap(Point1, Point0); }
    if (Point2.y < Point0.y) { std::swap(Point2, Point0); }
    if (Point2.y < Point1.y) { std::swap(Point2, Point1); }

    // Compute the x coordinates of the triangle edges
    std::vector<int> x01 = Interpolate(Point0, Point1);
    std::vector<int> x12 = Interpolate(Point1, Point2);
    std::vector<int> x02 = Interpolate(Point0, Point2);

    x01.pop_back();
    x12.pop_back();

    auto x012 = x01;
    x012.insert( x012.end(), x12.begin(), x12.end() );

    int MiddleIndex = (int)(x012.size()) / 2;
    std::vector<int> Left, Right;
    if (x012[MiddleIndex] < x02[MiddleIndex]) {
        Left = x012;
        Right = x02;
    } else {
        Left = x02;
        Right = x012;
    }

    int TopY = (int)(Point0.y);
    int Iterations = (int)(x012.size());
    for (int i = 0; i < Iterations; i++) {
        for (int x = Left[i]; x <= Right[i]; x++)
            SDL_RenderDrawPoint(Renderer, x, TopY + i);
    }
}

int main() {
    SDL_version Version;
    SDL_GetVersion(&Version);
    std::cout << "Using SDL " << (int)Version.major << "." << (int)Version.minor << "." << (int)Version.patch << std::endl;
    if (SDL_Init(SDL_INIT_EVERYTHING))
        std::cout << SDL_GetError() << std::endl;

    int NumOfDisplays = SDL_GetNumVideoDisplays();
    std::cout << "Found " << NumOfDisplays << " display/s" << std::endl;

    int LargestRefreshRate = -1;
    for (int i = 0; i < NumOfDisplays; i++) {
        SDL_DisplayMode DisplayInfo;
        SDL_GetDesktopDisplayMode(i, &DisplayInfo);
        int DisplayRefreshRate = DisplayInfo.refresh_rate;
        if (DisplayRefreshRate > LargestRefreshRate)
            LargestRefreshRate = DisplayRefreshRate;
        int DisplayWidth = DisplayInfo.w;
        int DisplayHeight = DisplayInfo.h;
        std::cout << " #" << i << ": " << DisplayWidth << "x" << DisplayHeight << " " << DisplayRefreshRate << "fps" << std::endl;
    }
    std::cout << "Using framerate " << LargestRefreshRate << std::endl;

    SDL_Window* Window = SDL_CreateWindow("CustomThreeDeeSDL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 600, 600, SDL_WINDOW_SHOWN);
    SDL_Renderer* Renderer = SDL_CreateRenderer(Window, -1, 0);
    Clock Timer( (double)(LargestRefreshRate) );

    int DisplayWidth = 128;
    int DisplayHeight = 128;
    SDL_Texture* Texture = SDL_CreateTexture(Renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_TARGET, DisplayWidth, DisplayHeight);

    SDL_SetRenderDrawBlendMode(Renderer, SDL_BLENDMODE_BLEND);

    mesh Cube;
    Cube.triangles = {
            // SOUTH
            { -0.5f, -0.5f, -0.5f,    -0.5f, 0.5f, -0.5f,    0.5f, 0.5f, -0.5f },
            { -0.5f, -0.5f, -0.5f,    0.5f, 0.5f, -0.5f,    0.5f, -0.5f, -0.5f },
            // EAST
            { 0.5f, -0.5f, -0.5f,    0.5f, 0.5f, -0.5f,    0.5f, 0.5f, 0.5f },
            { 0.5f, -0.5f, -0.5f,    0.5f, 0.5f, 0.5f,    0.5f, -0.5f, 0.5f },
            // NORTH
            { 0.5f, -0.5f, 0.5f,    0.5f, 0.5f, 0.5f,    -0.5f, 0.5f, 0.5f },
            { 0.5f, -0.5f, 0.5f,    -0.5f, 0.5f, 0.5f,    -0.5f, -0.5f, 0.5f },
            // WEST
            { -0.5f, -0.5f, 0.5f,    -0.5f, 0.5f, 0.5f,    -0.5f, 0.5f, -0.5f },
            { -0.5f, -0.5f, 0.5f,    -0.5f, 0.5f, -0.5f,    -0.5f, -0.5f, -0.5f },
            // TOP
            { -0.5f, 0.5f, -0.5f,    -0.5f, 0.5f, 0.5f,    0.5f, 0.5f, 0.5f },
            { -0.5f, 0.5f, -0.5f,    0.5f, 0.5f, 0.5f,    0.5f, 0.5f, -0.5f },
            // BOTTOM
            { 0.5f, -0.5f, 0.5f,    -0.5f, -0.5f, 0.5f,    -0.5f, -0.5f, -0.5f },
            { 0.5f, -0.5f, 0.5f,    -0.5f, -0.5f, -0.5f,    0.5f, -0.5f, -0.5f },
    };

    vec3 Camera = { 0, 0, 0 };
    vec3 LightDirection = { 0, 0, -1 };
    double LightLength = std::sqrt(std::pow(LightDirection.x, 2) + std::pow(LightDirection.y, 2) + std::pow(LightDirection.z, 2));
    LightDirection.x /= LightLength;
    LightDirection.y /= LightLength;
    LightDirection.z /= LightLength;

    Matrix ProjectionMatrix;
    SetProjectionMatrix(ProjectionMatrix, 0.1, 1000.0, 90.0, (double)DisplayWidth / (double)DisplayHeight);

    bool Running = true;
    while (Running) {
        SDL_Event Event;
        while (SDL_PollEvent(&Event)) {
            switch (Event.type) {
                case SDL_QUIT: {
                    Running = false;
                } break;
            }
        }

        double TimeElapsed = Timer.GetTotalTimeElapsed();

        SDL_SetRenderTarget(Renderer, Texture);
        SDL_SetRenderDrawColor(Renderer, 0, 0, 0, 255);
        SDL_RenderClear(Renderer);

        for (TrianglePlane &triangle : Cube.triangles) {
            Matrix RotationXMatrix;
            SetRotationXMatrix(RotationXMatrix, TimeElapsed * 1.6);
            auto RotatedXTrianglePlane = TrianglePlane();
            PointApplyMatrix(triangle.points[0], RotatedXTrianglePlane.points[0], RotationXMatrix);
            PointApplyMatrix(triangle.points[1], RotatedXTrianglePlane.points[1], RotationXMatrix);
            PointApplyMatrix(triangle.points[2], RotatedXTrianglePlane.points[2], RotationXMatrix);

            Matrix RotationZMatrix;
            SetRotationZMatrix(RotationZMatrix, TimeElapsed * 1.0);
            auto RotatedZTrianglePlane = TrianglePlane();
            PointApplyMatrix(RotatedXTrianglePlane.points[0], RotatedZTrianglePlane.points[0], RotationZMatrix);
            PointApplyMatrix(RotatedXTrianglePlane.points[1], RotatedZTrianglePlane.points[1], RotationZMatrix);
            PointApplyMatrix(RotatedXTrianglePlane.points[2], RotatedZTrianglePlane.points[2], RotationZMatrix);

            Matrix RotationYMatrix;
            SetRotationYMatrix(RotationYMatrix, TimeElapsed * 0.9);
            auto RotatedYTrianglePlane = TrianglePlane();
            PointApplyMatrix(RotatedZTrianglePlane.points[0], RotatedYTrianglePlane.points[0], RotationYMatrix);
            PointApplyMatrix(RotatedZTrianglePlane.points[1], RotatedYTrianglePlane.points[1], RotationYMatrix);
            PointApplyMatrix(RotatedZTrianglePlane.points[2], RotatedYTrianglePlane.points[2], RotationYMatrix);

            auto MovedTrianglePlane = TrianglePlane();
            MovedTrianglePlane = RotatedYTrianglePlane;
            MovedTrianglePlane.points[0].z += 1.5;
            MovedTrianglePlane.points[1].z += 1.5;
            MovedTrianglePlane.points[2].z += 1.5;

            vec3 normal = vec3(), start = vec3(), end = vec3();
            start.x = MovedTrianglePlane.points[1].x - MovedTrianglePlane.points[0].x;
            start.y = MovedTrianglePlane.points[1].y - MovedTrianglePlane.points[0].y;
            start.z = MovedTrianglePlane.points[1].z - MovedTrianglePlane.points[0].z;
            end.x = MovedTrianglePlane.points[2].x - MovedTrianglePlane.points[0].x;
            end.y = MovedTrianglePlane.points[2].y - MovedTrianglePlane.points[0].y;
            end.z = MovedTrianglePlane.points[2].z - MovedTrianglePlane.points[0].z;
            normal.x = start.y * end.z - start.z * end.y;
            normal.y = start.z * end.x - start.x * end.z;
            normal.z = start.x * end.y - start.y * end.x;

            double NormalLength = std::sqrt(std::pow(normal.x, 2) + std::pow(normal.y, 2) + std::pow(normal.z, 2));
            normal.x /= NormalLength;
            normal.y /= NormalLength;
            normal.z /= NormalLength;

            double MagicNumberIDKWhatItDoes = normal.x * (MovedTrianglePlane.points[0].x - Camera.x) +
                normal.y * (MovedTrianglePlane.points[0].y - Camera.y) +
                normal.z * (MovedTrianglePlane.points[0].z - Camera.z);
            if (MagicNumberIDKWhatItDoes <= 0.0) {
                double Light = normal.x * LightDirection.x +
                                normal.y * LightDirection.y +
                                normal.z * LightDirection.z;

                auto ProjectedTrianglePlane = TrianglePlane();
                PointApplyMatrix(MovedTrianglePlane.points[0], ProjectedTrianglePlane.points[0], ProjectionMatrix);
                PointApplyMatrix(MovedTrianglePlane.points[1], ProjectedTrianglePlane.points[1], ProjectionMatrix);
                PointApplyMatrix(MovedTrianglePlane.points[2], ProjectedTrianglePlane.points[2], ProjectionMatrix);

                Triangle ProjectedTriangle(ProjectedTrianglePlane);
                PointMove(ProjectedTriangle.m_aPoints[0], 1.0, 1.0);
                PointMove(ProjectedTriangle.m_aPoints[1], 1.0, 1.0);
                PointMove(ProjectedTriangle.m_aPoints[2], 1.0, 1.0);

                PointScale(ProjectedTriangle.m_aPoints[0], 0.5 * DisplayWidth, 0.5 * DisplayHeight);
                PointScale(ProjectedTriangle.m_aPoints[1], 0.5 * DisplayWidth, 0.5 * DisplayHeight);
                PointScale(ProjectedTriangle.m_aPoints[2], 0.5 * DisplayWidth, 0.5 * DisplayHeight);

                ColorRGB Color = HSVtoRGB({(double) SDL_GetTicks() / 10.0, 0.5, Light});
                SDL_SetRenderDrawColor(Renderer, Color.r, Color.g, Color.b, 255);
                FillTriangle(Renderer, ProjectedTriangle);

                Color = HSVtoRGB({(double) SDL_GetTicks() / 10.0, 1.0, Light});
                SDL_SetRenderDrawColor(Renderer, Color.r, Color.g, Color.b, 255);
                DrawTriangle(Renderer, ProjectedTriangle);
            }
        }

        SDL_SetRenderTarget(Renderer, nullptr);
        SDL_RenderCopy(Renderer, Texture, nullptr, nullptr);
        SDL_RenderPresent(Renderer);
        Timer.Tick();
    }

    SDL_DestroyTexture(Texture);
    SDL_DestroyRenderer(Renderer);
    SDL_DestroyWindow(Window);
    SDL_Quit();
    return 0;
}

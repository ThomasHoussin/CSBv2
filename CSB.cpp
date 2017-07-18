#include <iostream>
#include <cmath>
#include <vector>

#define PI 3.14159265358979323846

class Point {
public:
	Point() {
		x = 0. ;
		y = 0. ;
	}
	Point(float x, float y) {
		this->x = x ;
		this->y = y ;
	}
	Point(Point const &p) {
		this->x = p.x ;
		this->y = p.y ;
	}
	float distance(float x, float y) {
		return sqrt((this->x - x) * (this->x - x) + (this->y - y) * (this->y - y)) ;
	}
	float distance(Point &p) {
		return this->distance(p.x, p.y) ;
	}
	float distanceSq(float x, float y) {
		return (this->x - x) * (this->x - x) + (this->y - y) * (this->y - y) ;
	}
	float distanceSq(Point &p) {
		return this->distanceSq(p.x, p.y) ;
	}

	float getX() const {
		return x;
	}

	void setX(float x) {
		this->x = x;
	}

	float getY() const {
		return y;
	}

	void setY(float y) {
		this->y = y;
	}

protected:
	float x ;
	float y ;
};

class Unit:public Point {
public:
	Unit(float x, float y):
		Point(x, y)
	{
		r = 0 ;
		vx = 0 ;
		vy = 0 ;
	}

	Unit(float x, float y, float vx, float vy):
		Point(x, y)
	{
		r = 0 ;
		this->vx = vx ;
		this->vy = vy ;
	}

	float getR() const {
		return r;
	}

	float getVx() const {
		return vx;
	}

	void setVx(float vx) {
		this->vx = vx;
	}

	float getVy() const {
		return vy;
	}

	void setVy(float vy) {
		this->vy = vy;
	}

protected:
	float r ;
	float vx ;
	float vy ;
} ;

class Checkpoint: public Unit {
public:
	Checkpoint(float x, float y):
		Unit(x, y)
	{
		this->r = 600 ;
	}
};

class Pod: public Unit {
public:
	Pod(float x, float y, float vx, float vy, float angle, int nextCheckpointId):
		Unit(x,y,vx,vy)
	{
		this->r = 400 ;
		this->angle = angle ;
		this->nextCheckpointId = nextCheckpointId ;
		this->timeout = 100 ;
		this->shield = false ;
	}

	float getAngle() const {
		return angle;
	}

	void setAngle(float angle) {
		this->angle = angle;
	}

	int getNextCheckpointId() const {
		return nextCheckpointId;
	}

	void setNextCheckpointId(int nextCheckpointId) {
		this->nextCheckpointId = nextCheckpointId;
	}

protected:
	float angle ;
	int nextCheckpointId ;
	int timeout ;
	bool shield ;
};

class Move {
public:
	Move(float angle, int thrust){
		this->angle = angle ;
		this->thrust = thrust ;
	}
	Move(Move const &m) {
		this->angle = m.angle ;
		this->thrust = m.thrust ;
	}
	float getAngle() const {
		return angle;
	}

	void setAngle(float angle) {
		this->angle = angle;
	}

	int getThrust() const {
		return thrust;
	}

	void setThrust(int thrust) {
		this->thrust = thrust;
	}

private :
	float angle; // Between -18.0 and +18.0
    int thrust; // Between -1 and 200
};

class Game {
private:
	static int laps ;
	static int checkpointCount ;
	static std::vector<Checkpoint*> checkpoints ;
	std::vector<Pod*> pods ;

public:
	Game(int laps,int checkpointCount) {
	    for (int i = 0; i < checkpointCount; i++) {
	        Game::laps = laps ;
	        Game::checkpointCount = checkpointCount ;
	    	int checkpointX;
	        int checkpointY;
	        std::cin >> checkpointX >> checkpointY; std::cin.ignore();
	        Checkpoint* c = new Checkpoint(checkpointX, checkpointY) ;
	        checkpoints.push_back(c) ;
	    }
	}

	void updatePods(int i, float x, float y, float vx, float vy, float angle, int nextCheckPointId) {
		if(pods.size() < 4) {
            Pod* pod = new Pod(x,y,vx,vy,angle,nextCheckPointId) ;
            pods.push_back(pod);
		}
		else {
			Pod* pod = pods[i] ;
			pod->setX(x);
			pod->setY(y);
			pod->setVx(vx);
			pod->setVy(vy);
			pod->setAngle(angle);
			pod->setNextCheckpointId(nextCheckPointId);
		}
	}

	void readGame() {
        for (int i = 0; i < 4; i++) {
            int x; // x position of your pod
            int y; // y position of your pod
            int vx; // x speed of your pod
            int vy; // y speed of your pod
            int angle; // angle of your pod
            int nextCheckPointId; // next check point id of your pod
            std::cin >> x >> y >> vx >> vy >> angle >> nextCheckPointId; std::cin.ignore();
            updatePods(i, x, y, vx, vy, angle, nextCheckPointId);
        }
	}

	static const int getCheckpointCount() {
		return checkpointCount;
	}

	static const int getLaps() {
		return laps;
	}
};

//static variables initialization
int Game::laps = 0 ;
int Game::checkpointCount = 0 ;
std::vector<Checkpoint*> Game::checkpoints ;

int main()
{
    int laps;
    std::cin >> laps; std::cin.ignore();
    int checkpointCount;
    std::cin >> checkpointCount; std::cin.ignore();

    Game* game = new Game(laps,checkpointCount) ;

    // game loop
    while (1) {
    	game->readGame();

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;


        // You have to output the target position
        // followed by the power (0 <= power <= 200)
        // i.e.: "x y power"
        std::cout << "8000 4500 100" << std::endl;
        std::cout << "8000 4500 100" << std::endl;
    }
}

//based on http://files.magusgeek.com/csb/csb_en.html

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>

#define PI 3.14159265358979323846

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

class Point {
public:
	float x ;
	float y ;

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
	float distance(const Point &p) {
		return this->distance(p.x, p.y) ;
	}
	float distanceSq(float x, float y) {
		return (this->x - x) * (this->x - x) + (this->y - y) * (this->y - y) ;
	}
	float distanceSq(const Point &p) {
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

	//It allows us to find the closest point on a line (described by two points) to another point.
	Point closest(const Point &a, const Point &b) {
	    float da = b.y - a.y;
	    float db = a.x - b.x;
	    float c1 = da*a.x + db*a.y;
	    float c2 = -db*this->x + da*this->y;
	    float det = da*da + db*db;
	    float cx = 0;
	    float cy = 0;

	    if (det != 0) {
	        cx = (da*c1 - db*c2) / det;
	        cy = (da*c2 + db*c1) / det;
	    } else {
	        // The point is already on the line
	        cx = this->x;
	        cy = this->y;
	    }

	    return Point(cx, cy);
	}

protected:

};

//forward declarations
class Unit ;
class Pod ;

class Collision {
public:
	Unit* a ;
	Unit* b;
	float t ;

	Collision(Unit* a, Unit* b, float t) :
		a(a), b(b), t(t) {

	}

	int operator==(Collision &object){
		if(this == &object) return true;
	    return false ;
	}
	int operator!=(Collision &object){
		if(this == &object) return false;
	    return true ;
	}

	static Collision NO_COLLISION ;
private:

};

Collision Collision::NO_COLLISION(NULL,NULL,2.) ;

class Unit:public Point {
public:
	Unit(float x, float y):
		Point(x, y)
	{
		r = 0 ;
		vx = 0 ;
		vy = 0 ;
	}

	virtual ~Unit() {
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

	virtual void bounce(Unit& unit) {
		throw std::logic_error("Appel à la fonction virtuelle Bounce de la classe Unit") ;
	}

	Collision collision(Unit* u) {
		// Square of the distance
	    float dist = this->distanceSq(*u);

	    // Sum of the radii squared
	    float sr = (this->r + u->r)*(this->r + u->r);

	    // We take everything squared to avoid calling sqrt uselessly. It is better for performances

	    if (dist < sr) {
	        // Objects are already touching each other. We have an immediate collision.
	        return Collision(this, u, 0.0);
	    }

	    // Optimisation. Objects with the same speed will never collide
	    if (this->vx == u->vx && this->vy == u->vy) {
	        return Collision::NO_COLLISION;
	    }

	    // We place ourselves in the reference frame of u. u is therefore stationary and is at (0,0)
	    float x = this->x - u->x;
	    float y = this->y - u->y;
	    Point myp(x, y);
	    float vx = this->vx - u->vx;
	    float vy = this->vy - u->vy;
	    Point up(0., 0.) ;

	    // We look for the closest point to u (which is in (0,0)) on the line described by our speed vector
	    Point p = up.closest(myp, Point(x + vx, y + vy));

	    // Square of the distance between u and the closest point to u on the line described by our speed vector
	    float pdist = up.distanceSq(p);

	    // Square of the distance between us and that point
	    float mypdist = myp.distanceSq(p);

	    // If the distance between u and this line is less than the sum of the radii, there might be a collision
	    if (pdist < sr) {
	     // Our speed on the line
	        float length = sqrt(vx*vx + vy*vy);

	        // We move along the line to find the point of impact
	        float backdist = sqrt(sr - pdist);
	        p.x = p.x - backdist * (vx / length);
	        p.y = p.y - backdist * (vy / length);

	        // If the point is now further away it means we are not going the right way, therefore the collision won't happen
	        if (myp.distanceSq(p) > mypdist) {
	            return Collision::NO_COLLISION;
	        }

	        pdist = p.distance(myp);

	        // The point of impact is further than what we can travel in one turn
	        if (pdist > length) {
	            return Collision::NO_COLLISION;
	        }

	        // Time needed to reach the impact point
	        float t = pdist / length;

	        return Collision(this, u, t);
	    }

	    return Collision::NO_COLLISION;
	}


protected:
	float r ;
	float vx ;
	float vy ;
} ;



class Pod: public Unit {
public:
	Pod(float x, float y, float vx, float vy, float angle, int nextCheckpointId):
		Unit(x,y,vx,vy)
	{
		this->r = 400 ;
		this->angle = angle ;
		this->nextCheckpointId = nextCheckpointId ;
		this->timeout = 100 ;
		this->shield = 0 ;
		this->checked = 0 ;
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

	bool isShieldActive() {
		return this->shield == 3 ;
	}

	void checkedCP() {
		this->timeout = 100 ;
		this->checked ++ ;
		this->nextCheckpointId ++ ;
	}

	float computeAngle(const Point& p) {
	    float d = this->distance(p);
	    float dx = (p.x - this->x) / d;
	    float dy = (p.y - this->y) / d;

	    // Simple trigonometry. We multiply by 180.0 / PI to convert radiants to degrees.
	    float a = acos(dx) * 180.0 / PI;

	    // If the point I want is below me, I have to shift the angle for it to be correct
	    if (dy < 0) {
	        a = 360.0 - a;
	    }
	    return a;
	}

	float diffAngle(const Point& p) {
	    float a = this->computeAngle(p);

	    // To know whether we should turn clockwise or not we look at the two ways and keep the smallest
	    // The ternary operators replace the use of a modulo operator which would be slower
	    float right = this->angle <= a ? a - this->angle : 360.0 - this->angle + a;
	    float left = this->angle >= a ? this->angle - a : this->angle + 360.0 - a;

	    if (right < left) {
	        return right;
	    } else {
	        // We return a negative angle if we must rotate to left
	        return -left;
	    }
	}

	void rotate(const Point& p) {
	    float a = this->diffAngle(p);

	    // Can't turn by more than 18° in one turn
	    if (a > 18.0) {
	        a = 18.0;
	    } else if (a < -18.0) {
	        a = -18.0;
	    }

	    this->angle += a;

	    // The % operator is slow. If we can avoid it, it's better.
	    if (this->angle >= 360.0) {
	        this->angle = this->angle - 360.0;
	    } else if (this->angle < 0.0) {
	        this->angle += 360.0;
	    }
	}

	void rotate(float angle) {
		if(angle > 18.0) angle = 18.0 ;
		else if(angle < -18.0) angle = -18.0 ;

		this->angle += angle;

		// The % operator is slow. If we can avoid it, it's better.
		if (this->angle >= 360.0) {
		    this->angle = this->angle - 360.0;
		} else if (this->angle < 0.0) {
		    this->angle += 360.0;
		}

	}

	void boost(int thrust) {
	  // Don't forget that a pod which has activated its shield cannot accelerate for 3 turns
	    if (this->shield > 0) {
	        return;
	    }

	    // Conversion of the angle to radiants
	    float ra = this->angle * PI / 180.0;

	    // Trigonometry
	    this->vx += cos(ra) * thrust;
	    this->vy += sin(ra) * thrust;
	}

	void move(float t) {
	    this->x += this->vx * t;
	    this->y += this->vy * t;
	}

	void endTurn(int checkpointCount) {
	    this->x = round(this->x);
	    this->y = round(this->y);
	    this->vx = (int)(this->vx * 0.85);
	    this->vy = (int)(this->vy * 0.85);

	    // Don't forget that the timeout goes down by 1 each turn. It is reset to 100 when you pass a checkpoint
	    this->timeout -= 1;
	    if(this->shield > 0) this->shield -= 1 ;
		if(this->nextCheckpointId == checkpointCount) this->nextCheckpointId = 0 ;
	}

	void bounce(Unit& unit) {
		Pod& u = (Pod&)unit ;

		// If a pod has its shield active its mass is 10 otherwise it's 1
		float m1 = this->isShieldActive() ? 10 : 1;
		float m2 = u.isShieldActive() ? 10 : 1;
		float mcoeff = (m1 + m2) / (m1 * m2);

		float nx = this->x - u.x;
		float ny = this->y - u.y;

		// Square of the distance between the 2 pods. This value could be hardcoded because it is always 800²
		//TODO
		float nxnysquare = nx*nx + ny*ny;

		float dvx = this->vx - u.vx;
		float dvy = this->vy - u.vy;

		// fx and fy are the components of the impact vector. product is just there for optimisation purposes
		float product = nx*dvx + ny*dvy;
		float fx = (nx * product) / (nxnysquare * mcoeff);
		float fy = (ny * product) / (nxnysquare * mcoeff);

		// We apply the impact vector once
		this->vx -= fx / m1;
		this->vy -= fy / m1;
		u.vx += fx / m2;
		u.vy += fy / m2;

		// If the norm of the impact vector is less than 120, we normalize it to 120
		float impulse = sqrt(fx*fx + fy*fy);
		if (impulse < 120.0) {
			fx = fx * 120.0 / impulse;
			fy = fy * 120.0 / impulse;
		}

		// We apply the impact vector a second time
		this->vx -= fx / m1;
		this->vy -= fy / m1;
		u.vx += fx / m2;
		u.vy += fy / m2;
	}

	void Output(Move& move) {
	    std::string s ;
		float a = angle + move.getAngle();

	    if (a >= 360.0) {
	        a = a - 360.0;
	    } else if (a < 0.0) {
	        a += 360.0;
	    }

	    // Look for a point corresponding to the angle we want
	    // Multiply by 10000.0 to limit rounding errors
	    a = a * PI / 180.0;
	    float px = this->x + cos(a) * 10000.0;
	    float py = this->y + sin(a) * 10000.0;

	    if (move.getThrust() == -1) {
	        std::cout << round(px) << " " << round(py) << " SHIELD" << std::endl;
	    } else {
	        std::cout << round(px) << " " << round(py) <<" " << move.getThrust() << std::endl;
	    }
	}

protected:
	float angle ;
	int nextCheckpointId ;
	int timeout ;
	int shield ;
	int checked ;
};

class Checkpoint: public Unit {
public:
	Checkpoint(float x, float y):
		Unit(x, y)
	{
		this->r = 200 ; //600 - 400
	}

	void bounce(Unit& unit) {
		((Pod&)unit).checkedCP() ;
	}
};

class Game {
private:
	static int laps ;
	static int checkpointCount ;

public:
	int turn ;
	static std::vector<Checkpoint> checkpoints ;
	std::vector<Pod> pods ;

	Game(int laps,int checkpointCount) {
	    turn = 0 ;

		for (int i = 0; i < checkpointCount; i++) {
	        Game::laps = laps ;
	        Game::checkpointCount = checkpointCount ;
	    	int checkpointX;
	        int checkpointY;
	        std::cin >> checkpointX >> checkpointY; std::cin.ignore();
	        Checkpoint c = Checkpoint(checkpointX, checkpointY) ;
	        checkpoints.push_back(c) ;
	    }
	}

	Game(const Game& g) {
		turn = g.turn ;

		for(Pod p:g.pods) {
			Pod cp = p ;
			pods.push_back(cp) ;
		}
	}

	void updatePods(int i, float x, float y, float vx, float vy, float angle, int nextCheckPointId) {
		if(pods.size() < 4) {
            Pod pod = Pod(x,y,vx,vy,angle,nextCheckPointId) ;
            pods.push_back(pod);
		}
		else {
			Pod& pod = pods[i] ;
			pod.setX(x);
			pod.setY(y);
			pod.setVx(vx);
			pod.setVy(vy);
			pod.setAngle(angle);
			pod.setNextCheckpointId(nextCheckPointId);
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

	void simulateMovement() {
	    // This tracks the time during the turn. The goal is to reach 1.0
	    float t = 0.0;

	    while (t < 1.0) {
	        Collision& firstCollision = Collision::NO_COLLISION ;
	        std::cerr << &Collision::NO_COLLISION << std::endl ;
	        std::cerr << &firstCollision << std::endl ;

	        // We look for all the collisions that are going to occur during the turn
	        for (int i = 0; i < static_cast<int>(pods.size()); i++) {
	            // Collision with another pod?
	            for (int j = i + 1; j < static_cast<int>(pods.size()); j++) {
	                Collision col = pods[i].collision(&pods[j]);

	                // If the collision occurs earlier than the one we currently have we keep it
	                if (col != Collision::NO_COLLISION && col.t + t < 1.0 && col.t < firstCollision.t) {
	                    firstCollision = col;
	                }
	            }

	            // Collision with another checkpoint?
	            // It is unnecessary to check all checkpoints here. We only test the pod's next checkpoint.
	            // We could look for the collisions of the pod with all the checkpoints, but if such a collision happens it wouldn't impact the game in any way
	            Collision col = pods[i].collision(&checkpoints[pods[i].getNextCheckpointId()]);

	            // If the collision happens earlier than the current one we keep it
	            if (col != Collision::NO_COLLISION && col.t + t < 1.0 && col.t < firstCollision.t) {
	                firstCollision = col;
	            }
	        }

	        if (firstCollision == Collision::NO_COLLISION) {
	            // No collision, we can move the pods until the end of the turn
	            for (int i = 0; i < static_cast<int>(pods.size()); i++) {
	                pods[i].move(1.0 - t);
	            }

	            // End of the turn
	            t = 1.0;
	        } else {
	            // Move the pods to reach the time `t` of the collision
	            for (int i = 0; i < static_cast<int>(pods.size()); ++i) {
	                pods[i].move(firstCollision.t);
	            }

	            // Play out the collision
	            firstCollision.b->bounce(*firstCollision.a) ;

	            t += firstCollision.t;
	        }
	    }

	    for (int i = 0; i < static_cast<int>(pods.size()); ++i) {
	        pods[i].endTurn(Game::getCheckpointCount());
	    }
	}

	void simulateNextTurn(const std::vector<Move>& moves) {
		for(int i = 0; i < static_cast<int>(moves.size()) ; i++) {
			this->pods[i].rotate(moves[i].getAngle());
			this->pods[i].boost(moves[i].getThrust()) ;
		}
		this->simulateMovement();
	}
};

//static variables initialization
int Game::laps = 0 ;
int Game::checkpointCount = 0 ;
std::vector<Checkpoint> Game::checkpoints ;

class IA {
public:
	virtual std::vector<Move> computeMoves() = 0 ;
	virtual ~IA() { } ;

	Game* game ;
	clock_t begin ;

	void resetTimer() {
		begin = clock() ;
	}

	bool timeOut() {
		int time = game->turn == 1 ? MAX_TIME_FIRST_TURN : MAX_TIME ;
		return  (float)(clock()-begin)/CLOCKS_PER_SEC * 1000 >= time ;
	}
private :
	static const int MAX_TIME_FIRST_TURN = 1000 ;
	static const int MAX_TIME = 150 ;
};

class SimpleIA : public IA {
public:
	SimpleIA(Game* game) {
		this->game = game ;
		begin = clock() ;
	}

	~SimpleIA() {

	}

	std::vector<Move> computeMoves() {
		std::vector<Move> moves ;
		for(int i = 0 ; i <2 ; i++) {
			Checkpoint& CP = game->checkpoints[game->pods[i].getNextCheckpointId()] ;
			Move m(game->pods[i].diffAngle(CP) , 100) ;
			if(m.getAngle() > 45) m.setThrust(0) ;
			moves.push_back(m) ;
		}
		return moves ;
	}
private:

};

int main()
{
    int laps;
    std::cin >> laps; std::cin.ignore();
    int checkpointCount;
    std::cin >> checkpointCount; std::cin.ignore();

    Game game(laps,checkpointCount) ;
    SimpleIA ia(&game) ;

    // game loop
    while (1) {
    	ia.resetTimer();
    	game.readGame();

    	std::vector<Move> moves = ia.computeMoves() ;

    	Game clone(game) ;
    	clone.simulateNextTurn(moves);
    	std::cerr << "x : " << clone.pods[0].x << " y : " << clone.pods[0].y << std::endl ;
    	std::cerr << "vx : " << clone.pods[0].getVx() << " vy : " << clone.pods[0].getVy() << std::endl ;

    	std::cerr << "x : " << clone.pods[1].x << " y : " << clone.pods[1].y << std::endl ;
    	std::cerr << "vx : " << clone.pods[1].getVx() << " vy : " << clone.pods[1].getVy() << std::endl ;

        // Write an action using cout. DON'T FORGET THE "<< endl"
        // To debug: cerr << "Debug messages..." << endl;


        // You have to output the target position
        // followed by the power (0 <= power <= 200)
        // i.e.: "x y power"
    	game.pods[0].Output(moves[0]);
    	game.pods[1].Output(moves[1]);
    }
}

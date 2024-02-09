import pygame
import math
import constants


class Body(pygame.sprite.Sprite):
    def __init__(self, x, y, mass, radius, color):
        self.x = x
        self.y = y

        # mass in kg
        self.m = mass
        self.color = color
        self.r = radius
        self.x_vel = 0
        self.y_vel = 0
        self.collision = False

    def angle(self, other):

        # find the angle between two bodies

        dotProd = self.x * other.x + self.y * other.y

        vectmagn = math.sqrt(self.x ** 2 + self.y ** 2) * math.sqrt(other.x ** 2 + other.y ** 2)

        CosAngle = dotProd / vectmagn

        angleDeg = math.degrees(math.acos(CosAngle))
        # print(angleDeg)
        return angleDeg

    def force(self, other):

        # x and y of other body
        x2, y2 = other.x, other.y

        # distance between this body and other body
        distancex = x2 - self.x
        distancey = y2 - self.y
        distance = math.sqrt(distancex ** 2 + distancey ** 2)

        # force vectors on body
        force = (constants.G * self.m * other.m) / distance ** 2
        theta = math.atan2(distancey, distancex)
        forcex = math.cos(theta) * force
        forcey = math.sin(theta) * force

        return forcex, forcey

    def draw(self, win):
        x = self.x * constants.SCALE + constants.SCREEN_WIDTH / 2
        y = self.y * constants.SCALE + constants.SCREEN_HEIGHT / 2
        pygame.draw.circle(win, self.color, (x, y), self.r)

    def update(self, bodies):

        if not self.collision:
            totalfx = totalfy = 0

            for body in bodies:

                if self == body:
                    continue

                fx, fy = self.force(body)
                totalfx += fx
                totalfy += fy

            self.x_vel += totalfx / self.m * constants.TIME_STEP
            self.y_vel += totalfy / self.m * constants.TIME_STEP

            self.x += self.x_vel * constants.TIME_STEP
            self.y += self.y_vel * constants.TIME_STEP

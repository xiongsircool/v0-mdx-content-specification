"use client"

import { motion } from 'framer-motion'
import { useInView } from 'react-intersection-observer'
import {
  Users,
  BookOpen,
  Zap,
  Code,
  Trophy,
  Heart,
  Rocket,
  Sparkles
} from 'lucide-react'
import { AnimatedCard, ParallaxSection } from './animations'

interface FeaturesGridProps {
  siteStats?: {
    activeMembers: number
    techShares: number
    learnResources: number
    totalActivities: number
  }
}

const features = [
  {
    icon: Users,
    title: "自由交流",
    description: "没有门槛，没有压力，只有对生物信息的热爱和分享",
    color: "from-blue-500 to-purple-500",
    bgGradient: "from-blue-50 to-purple-50",
    hoverGradient: "from-blue-100 to-purple-100"
  },
  {
    icon: BookOpen,
    title: "学习分享",
    description: "一起学习新技术，分享学习资料，互助成长",
    color: "from-purple-500 to-pink-500",
    bgGradient: "from-purple-50 to-pink-50",
    hoverGradient: "from-purple-100 to-pink-100"
  },
  {
    icon: Code,
    title: "技术实践",
    description: "动手实践真实项目，在实践中提升技能",
    color: "from-pink-500 to-yellow-500",
    bgGradient: "from-pink-50 to-yellow-50",
    hoverGradient: "from-pink-100 to-yellow-100"
  },
  {
    icon: Rocket,
    title: "创新探索",
    description: "勇于尝试新技术，探索生物信息学的无限可能",
    color: "from-yellow-500 to-green-500",
    bgGradient: "from-yellow-50 to-green-50",
    hoverGradient: "from-yellow-100 to-green-100"
  },
  {
    icon: Heart,
    title: "社区温暖",
    description: "温暖的社区氛围，让每个人都能找到归属感",
    color: "from-green-500 to-blue-500",
    bgGradient: "from-green-50 to-blue-50",
    hoverGradient: "from-green-100 to-blue-100"
  },
  {
    icon: Trophy,
    title: "共同成长",
    description: "一起挑战，一起进步，见证彼此的成长",
    color: "from-blue-500 to-indigo-500",
    bgGradient: "from-blue-50 to-indigo-50",
    hoverGradient: "from-blue-100 to-indigo-100"
  }
]

export function FeaturesGrid({ siteStats }: FeaturesGridProps) {
  const [ref, inView] = useInView({
    triggerOnce: true,
    threshold: 0.1,
  })

  return (
    <section className="py-20 px-4 bg-gradient-to-br from-gray-50 to-white">
      <ParallaxSection speed={0.3}>
        <div className="container mx-auto">
          <div className="text-center mb-16">
            <motion.div
              ref={ref}
              initial={{ opacity: 0, y: 50 }}
              animate={inView ? { opacity: 1, y: 0 } : { opacity: 0, y: 50 }}
              transition={{ duration: 0.8 }}
            >
              <h2 className="text-4xl md:text-5xl font-bold mb-4 bg-gradient-to-r from-blue-600 via-purple-600 to-pink-600 bg-clip-text text-transparent">
                我们在这里做什么
              </h2>
              <p className="text-xl text-gray-600 max-w-3xl mx-auto">
                一群热爱生物信息的学生，一起学习、分享、搞事情
              </p>
            </motion.div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-8">
            {features.map((feature, index) => {
              const Icon = feature.icon
              return (
                <AnimatedCard key={feature.title} delay={index * 0.1}>
                  <motion.div
                    className={`group relative bg-gradient-to-br ${feature.bgGradient} rounded-3xl p-8 shadow-lg hover:shadow-2xl transition-all duration-500 overflow-hidden`}
                    whileHover={{
                      y: -10,
                      background: `linear-gradient(135deg, ${feature.hoverGradient})`,
                    }}
                  >
                    {/* Background decoration */}
                    <div className="absolute top-0 right-0 w-32 h-32 bg-gradient-to-br opacity-10 rounded-full blur-2xl" />

                    {/* Icon */}
                    <div className="flex items-center justify-center w-16 h-16 bg-gradient-to-r rounded-2xl mb-6 group-hover:scale-110 transition-transform duration-300">
                      <Icon className="w-8 h-8 text-white" />
                    </div>

                    {/* Content */}
                    <h3 className="text-2xl font-bold mb-4 text-gray-900 group-hover:text-transparent group-hover:bg-clip-text group-hover:bg-gradient-to-r group-hover:from-blue-600 group-hover:to-purple-600 transition-all duration-300">
                      {feature.title}
                    </h3>

                    <p className="text-gray-600 leading-relaxed group-hover:text-gray-700 transition-colors duration-300">
                      {feature.description}
                    </p>

                    {/* Floating sparkles */}
                    <motion.div
                      className="absolute bottom-4 right-4 opacity-0 group-hover:opacity-100 transition-opacity duration-300"
                      animate={{ rotate: 360 }}
                      transition={{ duration: 3, repeat: Infinity, ease: "linear" }}
                    >
                      <Sparkles className="w-5 h-5 text-blue-500" />
                    </motion.div>
                  </motion.div>
                </AnimatedCard>
              )
            })}
          </div>

          {/* Stats section */}
          <div className="mt-20 grid grid-cols-2 md:grid-cols-4 gap-8">
            {[
              {
                label: "活跃爱好者",
                value: siteStats ? `${siteStats.activeMembers}+` : "50+",
                color: "from-blue-500 to-purple-500"
              },
              {
                label: "技术分享",
                value: siteStats ? `${siteStats.techShares}+` : "20+",
                color: "from-purple-500 to-pink-500"
              },
              {
                label: "学习资料",
                value: siteStats ? `${siteStats.learnResources}+` : "10+",
                color: "from-pink-500 to-yellow-500"
              },
              {
                label: "社区活动",
                value: siteStats ? `${siteStats.totalActivities}+` : "10+",
                color: "from-yellow-500 to-green-500"
              },
            ].map((stat, index) => (
              <AnimatedCard key={stat.label} delay={index * 0.1 + 0.5}>
                <div className="text-center">
                  <div className={`text-4xl md:text-6xl font-bold bg-gradient-to-r ${stat.color} bg-clip-text text-transparent mb-2`}>
                    {stat.value}
                  </div>
                  <p className="text-lg text-gray-600 font-medium">{stat.label}</p>
                </div>
              </AnimatedCard>
            ))}
          </div>
        </div>
      </ParallaxSection>
    </section>
  )
}
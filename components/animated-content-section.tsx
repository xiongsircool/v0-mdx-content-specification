"use client"

import { motion } from 'framer-motion'
import { useInView } from 'react-intersection-observer'
import { ContentCard } from './content-card'
import type { ContentItem } from '@/lib/server-markdown-loader'
import { Megaphone, BookOpen, ChevronRight } from 'lucide-react'

interface AnimatedContentSectionProps {
  announcements: ContentItem[]
  posts: ContentItem[]
}

export function AnimatedContentSection({ announcements, posts }: AnimatedContentSectionProps) {
  const [announcementsRef, announcementsInView] = useInView({
    triggerOnce: true,
    threshold: 0.1,
  })

  const [postsRef, postsInView] = useInView({
    triggerOnce: true,
    threshold: 0.1,
  })

  return (
    <>
      {/* Recent Announcements */}
      <section className="py-16 px-4 bg-gradient-to-br from-gray-50 to-white">
        <div className="container mx-auto">
          <div className="flex items-center justify-between mb-8">
            <div className="flex items-center gap-3">
              <div className="flex h-12 w-12 items-center justify-center rounded-xl bg-gradient-to-br from-blue-600 to-purple-600">
                <Megaphone className="h-6 w-6 text-white" />
              </div>
              <div>
                <h2 className="text-3xl font-bold bg-gradient-to-r from-blue-600 to-purple-600 bg-clip-text text-transparent">活动公告</h2>
                <p className="text-muted-foreground">了解联盟最新活动和招新信息</p>
              </div>
            </div>
            <a
              href="/announcements"
              className="inline-flex items-center gap-2 px-4 py-2 bg-transparent border border-blue-200 text-blue-600 hover:bg-blue-50 rounded-lg transition-all duration-300"
            >
              查看全部 <ChevronRight className="w-4 h-4" />
            </a>
          </div>

          {announcements.length > 0 ? (
            <div ref={announcementsRef} className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {announcements.map((announcement, index) => (
                <motion.div
                  key={announcement.slug}
                  initial={{ opacity: 0, y: 30 }}
                  animate={announcementsInView ? { opacity: 1, y: 0 } : { opacity: 0, y: 30 }}
                  transition={{ duration: 0.5, delay: index * 0.1 }}
                >
                  <div className="group relative bg-white rounded-2xl shadow-lg hover:shadow-2xl transition-all duration-300 overflow-hidden">
                    <div className="absolute inset-0 bg-gradient-to-br from-blue-600/5 to-purple-600/5 opacity-0 group-hover:opacity-100 transition-opacity duration-300" />
                    <div className="relative p-6">
                      <ContentCard item={announcement} type="announcement" />
                    </div>
                  </div>
                </motion.div>
              ))}
            </div>
          ) : (
            <div className="bg-gradient-to-br from-gray-50 to-white border-0 shadow-lg rounded-2xl p-12 text-center">
              <Megaphone className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">暂无公告</h3>
              <p className="text-muted-foreground">目前没有发布任何公告</p>
            </div>
          )}
        </div>
      </section>

      {/* Recent Posts */}
      <section className="py-16 px-4 bg-gradient-to-br from-purple-50 to-pink-50">
        <div className="container mx-auto">
          <div className="flex items-center justify-between mb-8">
            <div className="flex items-center gap-3">
              <div className="flex h-12 w-12 items-center justify-center rounded-xl bg-gradient-to-br from-purple-600 to-pink-600">
                <BookOpen className="h-6 w-6 text-white" />
              </div>
              <div>
                <h2 className="text-3xl font-bold bg-gradient-to-r from-purple-600 to-pink-600 bg-clip-text text-transparent">技术分享</h2>
                <p className="text-muted-foreground">成员技术分享与项目成果</p>
              </div>
            </div>
            <a
              href="/posts"
              className="inline-flex items-center gap-2 px-4 py-2 bg-transparent border border-purple-200 text-purple-600 hover:bg-purple-50 rounded-lg transition-all duration-300"
            >
              查看全部 <ChevronRight className="w-4 h-4" />
            </a>
          </div>

          {posts.length > 0 ? (
            <div ref={postsRef} className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6">
              {posts.map((post, index) => (
                <motion.div
                  key={post.slug}
                  initial={{ opacity: 0, y: 30 }}
                  animate={postsInView ? { opacity: 1, y: 0 } : { opacity: 0, y: 30 }}
                  transition={{ duration: 0.5, delay: index * 0.1 }}
                >
                  <div className="group relative bg-white rounded-2xl shadow-lg hover:shadow-2xl transition-all duration-300 overflow-hidden">
                    <div className="absolute inset-0 bg-gradient-to-br from-purple-600/5 to-pink-600/5 opacity-0 group-hover:opacity-100 transition-opacity duration-300" />
                    <div className="relative p-6">
                      <ContentCard item={post} type="post" />
                    </div>
                  </div>
                </motion.div>
              ))}
            </div>
          ) : (
            <div className="bg-gradient-to-br from-purple-50 to-pink-50 border-0 shadow-lg rounded-2xl p-12 text-center">
              <BookOpen className="h-16 w-16 text-muted-foreground mx-auto mb-4" />
              <h3 className="text-xl font-semibold mb-2">暂无文章</h3>
              <p className="text-muted-foreground">目前没有发布任何技术文章</p>
            </div>
          )}
        </div>
      </section>
    </>
  )
}